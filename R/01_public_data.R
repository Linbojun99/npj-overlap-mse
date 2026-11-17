# 01_public_data.R — 读入 WCPFC 与 IATTC 公开延绳钓数据，整理成同一张长表
# 01_public_data.R — read the WCPFC and IATTC public longline data into one long table
#
# 输入 | inputs (下载后放在 data/raw/，见 README | download into data/raw/, see README):
#   WCPFC_L_PUBLIC_BY_YY_MM_FLAG.csv            WCPFC 年-月-船旗-5°格 | year-month-flag-5° cell
#   PublicLLTunaBillfishNum.csv / ...Mt.csv     IATTC 尾数与重量两个文件 | IATTC numbers and weights
# 输出 | output: output/01_ll_long.csv

source("R/00_settings.R")

# ---- 1. 检查原始文件 | check the raw files ----
files <- file.path(RAW, c("WCPFC_L_PUBLIC_BY_YY_MM_FLAG.csv",
                                   "PublicLLTunaBillfishNum.csv", "PublicLLTunaBillfishMt.csv"))
if (!all(file.exists(files)))
  stop("缺少原始数据，见 README | raw data missing; see README")

# ---- 2. 分析范围 | analysis domain ----
year_from <- 1990; year_to <- 2024
# 只保留整格在赤道以南的格（西南角纬度 ≤ −5）| keep cells entirely south of the equator (SW-corner latitude ≤ −5)
sp_lat_sw_max <- -5
# 重疊区：4°S–50°S，150°W–130°W（IATTC C-12-11）| overlap area: 4°S–50°S, 150°W–130°W (IATTC C-12-11)
ov_south <- -50; ov_north <- -4; ov_west <- -150; ov_east <- -130

# ---- 3. 船旗对照 | flag crosswalk (WCPFC ISO2 ↔ IATTC ISO3) ----
flag_xw <- tribble(
  ~iso2, ~iso3, ~flag_label,
  "JP", "JPN", "Japan",            "TW", "TWN", "Chinese Taipei",  "KR", "KOR", "Korea",
  "CN", "CHN", "China",            "US", "USA", "USA",             "VU", "VUT", "Vanuatu",
  "PF", "PYF", "French Polynesia", "FJ", "FJI", "Fiji",            "ES", "ESP", "Spain",
  "PA", "PAN", "Panama",           "BZ", "BLZ", "Belize",          "MX", "MEX", "Mexico",
  "NZ", "NZL", "New Zealand",      "AU", "AUS", "Australia",       "NC", "NCL", "New Caledonia",
  "CK", "COK", "Cook Islands",     "WS", "WSM", "Samoa",           "TO", "TON", "Tonga",
  "SB", "SLB", "Solomon Is.",      "FM", "FSM", "Micronesia",      "KI", "KIR", "Kiribati")

# ---- 4. WCPFC：西南角标签 "05S" "150W"，钩数以百钩计 ----
# ---- 4. WCPFC: SW-corner labels "05S" "150W"; hooks in hundreds ----
wcpfc_raw <- read_csv(files[1], col_types = cols(.default = col_character()))
lat_v <- as.numeric(substr(wcpfc_raw$LAT5, 1, nchar(wcpfc_raw$LAT5) - 1))   # 数字部分 | number part
lon_v <- as.numeric(substr(wcpfc_raw$LON5, 1, nchar(wcpfc_raw$LON5) - 1))
lat_h <- substr(wcpfc_raw$LAT5, nchar(wcpfc_raw$LAT5), nchar(wcpfc_raw$LAT5)) # 半球字母 | hemisphere letter
lon_h <- substr(wcpfc_raw$LON5, nchar(wcpfc_raw$LON5), nchar(wcpfc_raw$LON5))
wcpfc <- tibble(
  rfmo   = "WCPFC",
  year   = as.integer(wcpfc_raw$YY), month = as.integer(wcpfc_raw$MM),
  code   = wcpfc_raw$flag_code,
  lat_sw = ifelse(lat_h == "S", -lat_v, lat_v),           # 南纬为负 | south negative
  lon_sw = ifelse(lon_h == "W", -lon_v, lon_v),           # 西经为负 | west negative
  hooks  = as.numeric(wcpfc_raw$HHOOKS) * 100,            # 百钩 → 钩 | hundreds → hooks
  alb_n  = as.numeric(wcpfc_raw$ALB_N),                   # 长鳍鲔尾数 | albacore numbers
  alb_mt = as.numeric(wcpfc_raw$ALB_C))                   # 长鳍鲔吨 | albacore tonnes
wcpfc$flag_label <- flag_xw$flag_label[match(wcpfc$code, flag_xw$iso2)]

# ---- 5. IATTC：格中心坐标，尾数与重量两个文件按键一对一合并 ----
# ---- 5. IATTC: cell-centre coordinates; number and weight files joined one-to-one on their keys ----
iattc_n  <- read_csv(files[2], show_col_types = FALSE)
iattc_mt <- read_csv(files[3], show_col_types = FALSE)
keys <- c("Year", "Month", "Flag", "LatC5", "LonC5")
stopifnot(!anyDuplicated(iattc_n[keys]), !anyDuplicated(iattc_mt[keys]),
          nrow(anti_join(iattc_n, iattc_mt, by = keys)) == 0)        # 两文件同一组键 | same keys in both
joined <- left_join(select(iattc_n, all_of(keys), Hooks, ALBn),
                    select(iattc_mt, all_of(keys), Hooks_mt = Hooks, ALBmt), by = keys)
stopifnot(all(joined$Hooks == joined$Hooks_mt))                      # 钩数只计一次 | hooks enter once
iattc <- tibble(
  rfmo   = "IATTC",
  year   = as.integer(joined$Year), month = as.integer(joined$Month),
  code   = joined$Flag,
  lat_sw = joined$LatC5 - 2.5, lon_sw = joined$LonC5 - 2.5,          # 中心 → 西南角 | centre → SW corner
  hooks  = as.numeric(joined$Hooks),
  alb_n  = as.numeric(joined$ALBn),
  alb_mt = as.numeric(joined$ALBmt))
iattc$flag_label <- flag_xw$flag_label[match(iattc$code, flag_xw$iso3)]

# ---- 6. 合并、筛选年份与南半球 | combine, filter years and the southern domain ----
ll <- bind_rows(wcpfc, iattc) |>
  filter(year >= year_from, year <= year_to, lat_sw <= sp_lat_sw_max)
# 空白船旗 = 保密遮蔽；未在对照表中的保留代码 | blank flag = confidentiality mask; unmapped codes kept
ll$flag_label <- ifelse(is.na(ll$code) | trimws(ll$code) == "", "Flag masked",
                        ifelse(is.na(ll$flag_label), paste0("Unmapped: ", ll$code), ll$flag_label))

# ---- 7. 每格落在重疊区内的比例（三种边界处理）| fraction of each cell inside the overlap (three boundary treatments) ----
# 北界 4°S 切过 5°S–0° 那一排格；inner = 只算整格在内，intersecting = 相交即全算，
# area = 按球面面积比例（假设格内均匀分布）。正文用 area。
# The 4°S boundary cuts the 5°S–0° row; inner = whole cells only, intersecting = any overlap counts fully,
# area = spherical-area fraction (uniform within the cell). The Article uses area.
south <- pmax(ll$lat_sw, ov_south); north <- pmin(ll$lat_sw + 5, ov_north)
west  <- pmax(ll$lon_sw, ov_west);  east  <- pmin(ll$lon_sw + 5, ov_east)
intersects <- north > south & east > west
ll$w_inner        <- as.numeric(ll$lat_sw >= ov_south & ll$lat_sw + 5 <= ov_north &
                                ll$lon_sw >= ov_west  & ll$lon_sw + 5 <= ov_east)
ll$w_intersecting <- as.numeric(intersects)
ll$w_area <- ifelse(intersects, (east - west) / 5 *
  (sin(north * pi / 180) - sin(south * pi / 180)) /
  (sin((ll$lat_sw + 5) * pi / 180) - sin(ll$lat_sw * pi / 180)), 0)
ll$in_overlap <- ll$w_intersecting == 1
stopifnot(all(is.finite(as.matrix(ll[c("hooks", "alb_n", "alb_mt")]))))

# ---- 8. 看一下结果 | look at the result ----
print(count(ll, rfmo, name = "rows"))
print(ll |> filter(in_overlap) |> count(rfmo, name = "overlap_rows"))
cat("北界格的面积权重 | area weight of a northern boundary cell:",
    round(unique(ll$w_area[ll$lat_sw == -5 & ll$lon_sw == -150]), 6), "\n")
write_csv(ll, file.path(OUT, "01_ll_long.csv"))
cat("01 done: output/01_ll_long.csv\n")
