# 02_overlap_share.R — 重疊区份额、两组织记录的船旗比值与记录匹配（图 1、补充表 1–3）
# 02_overlap_share.R — overlap-area shares, flag ratios between the two products and record matching
#                      (Fig. 1, Supplementary Tables 1–3)
#
# 输入 | input: output/01_ll_long.csv（01 生成 | from step 01）
# 输出 | outputs: output/02_*.csv

source("R/00_settings.R")
if (!exists("ll")) ll <- read_csv(file.path(OUT, "01_ll_long.csv"), show_col_types = FALSE)

# ---- 1. 每年重疊区占各自产品赤道以南总量的份额 | annual overlap share of each product's south-of-equator total ----
# 分母是该组织自己的公开产品，不是全种群渔获 | denominator is each commission's own public product, not stock-wide catch
share_year <- bind_rows(
  ll |> mutate(scheme = "inner",        w = w_inner),
  ll |> mutate(scheme = "area",         w = w_area),
  ll |> mutate(scheme = "intersecting", w = w_intersecting)) |>
  group_by(rfmo, scheme, year) |>
  summarise(hooks_total = sum(hooks),  hooks_overlap = sum(hooks * w),
            albn_total  = sum(alb_n),  albn_overlap  = sum(alb_n * w),
            albmt_total = sum(alb_mt), albmt_overlap = sum(alb_mt * w), .groups = "drop") |>
  mutate(share_hooks = hooks_overlap / hooks_total,
         share_albn  = albn_overlap / albn_total,
         share_albmt = albmt_overlap / albmt_total)
write_csv(share_year, file.path(OUT, "02_overlap_share_by_year.csv"))

# ---- 2. 2015–2024 平均份额（补充表 1）| 2015–2024 shares (Supplementary Table 1) ----
period <- share_year |> filter(year >= 2015, year <= 2024) |> group_by(rfmo, scheme) |>
  summarise(n_years = n(),
            mean_annual_hooks_share = mean(share_hooks), mean_annual_mt_share = mean(share_albmt),
            pooled_hooks_share = sum(hooks_overlap) / sum(hooks_total),
            pooled_mt_share    = sum(albmt_overlap) / sum(albmt_total), .groups = "drop")
print(period, width = Inf)
write_csv(period, file.path(OUT, "02_period_shares.csv"))

# ---- 3. 按船旗：两组织 2015–2024 重疊区钩数、重量与尾数（图 1d、补充表 2）----
# ---- 3. by flag: 2015–2024 overlap hooks, weight and numbers in the two products (Fig. 1d, Supplementary Table 2) ----
flag_year <- ll |> filter(in_overlap, year >= 2015) |>
  group_by(rfmo, flag_label, year) |>
  summarise(hooks = sum(hooks * w_area), alb_mt = sum(alb_mt * w_area), alb_n = sum(alb_n * w_area), .groups = "drop")
flag_wide <- flag_year |>
  pivot_longer(c(hooks, alb_mt, alb_n), names_to = "metric") |>
  pivot_wider(names_from = rfmo, values_from = value, values_fill = 0)
flag_ratio <- flag_wide |> group_by(flag_label, metric) |>
  # 先数两个产品都为正值的年数，再求和（summarise 按顺序计算，求和会覆盖逐年列）
  # count years positive in both products first, then sum (summarise evaluates in order; the sums overwrite the annual columns)
  summarise(n_common_positive = sum(WCPFC > 0 & IATTC > 0),
            WCPFC = sum(WCPFC), IATTC = sum(IATTC), .groups = "drop") |>
  mutate(ratio = ifelse(WCPFC > 0 & IATTC > 0, IATTC / WCPFC, NA_real_),          # IATTC / WCPFC
         all_ten_years = n_common_positive == 10)
print(flag_ratio |> filter(metric == "hooks") |> arrange(desc(WCPFC)), width = Inf)
write_csv(flag_ratio, file.path(OUT, "02_flag_ratio_2015_2024.csv"))

# 正文：十年都有正值的五个船旗，IATTC 的十年钩数总量比 WCPFC 高 10–77%
# Article: for the five flags positive in both products in every year, IATTC decade hook totals exceed WCPFC by 10–77%
five <- flag_ratio |> filter(metric == "hooks", all_ten_years)
cat("五个船旗 | five flags:", paste(five$flag_label, collapse = ", "), "\n")
cat("钩数差异 | hook difference:", sprintf("%.0f–%.0f%%", 100 * (min(five$ratio) - 1),
                                           100 * (max(five$ratio) - 1)), "\n")
twn <- flag_ratio |> filter(metric == "alb_mt", flag_label == "Chinese Taipei")
cat("Chinese Taipei 长鳍鲔重量 IATTC/WCPFC | albacore weight IATTC/WCPFC:", round(twn$ratio, 3), "\n")

# ---- 4. 记录匹配：年-月-船旗-格，两产品都有正钩数的记录（正文 1,167 / 4,698 / 42）----
# ---- 4. record matching on year-month-flag-cell with positive hooks (Article: 1,167 / 4,698 / 42) ----
key <- c("year", "month", "flag_label", "lat_sw", "lon_sw")
w_keys <- ll |> filter(rfmo == "WCPFC", year >= 2015, in_overlap, hooks > 0) |> select(all_of(key))
i_keys <- ll |> filter(rfmo == "IATTC", year >= 2015, in_overlap, hooks > 0) |> select(all_of(key))
both <- nrow(inner_join(w_keys, i_keys, by = key))
support <- tibble(support = c("Both positive hooks", "IATTC only", "WCPFC only"),
                  n_keys = c(both, nrow(i_keys) - both, nrow(w_keys) - both))
print(support)
write_csv(support, file.path(OUT, "02_key_support.csv"))

# ---- 5. WCPFC 保密规则删去的钩数比例（整个公开文件）| hooks removed by WCPFC confidentiality (whole release) ----
cov <- read_csv(file.path(RAW, "WCPFC_L_PUBLIC_BY_YY_MM_FLAG_Coverage.csv"), show_col_types = FALSE) |>
  filter(yy >= 2015, yy <= 2024) |> mutate(removed_percent = 100 - Hks_reported_perc)
cat("各年删去 | removed by year:", sprintf("%.1f–%.1f%%", min(cov$removed_percent), max(cov$removed_percent)), "\n")

cat("02 done\n")
