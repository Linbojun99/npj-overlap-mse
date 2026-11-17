# 12_figures_main.R — 正文图 1–6 | Article Figs 1–6
# 需要 01–02、07（12 个配置）与 08 的结果 | needs the results of steps 01–02, 07 (12 configurations) and 08
# EEZ 外框为可选项：设定 EEZ_SHP 指向 Marine Regions EEZ v12（0–360 经度版）并安装 sf
# EEZ outlines are optional: set EEZ_SHP to the Marine Regions EEZ v12 shapefile (0–360 longitudes) and install sf

source("R/00_settings.R")
COL  <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7")      # MP-A 蓝、MP-B 橙、MP-C 绿、MP-D 粉 | blue, orange, green, pink
PCH  <- c(16, 15, 17, 18)                                  # 圆、方、三角、菱形 | circle, square, triangle, diamond
MPS  <- c("MP-A WCPFC TAC", "MP-B joint TAC", "MP-C TAC with deduction", "MP-D overlap effort cap")
LAB  <- c("MP-A", "MP-B", "MP-C", "MP-D")
A_DIR <- file.path(OUT, "08_analysis")
marg   <- read.csv(file.path(A_DIR, "marginal.csv"))
absd   <- read.csv(file.path(A_DIR, "absolute.csv"))
paired <- read.csv(file.path(A_DIR, "paired.csv"))
traj   <- read.csv(file.path(OUT, "07_reference", "first_replicate_trajectories.csv"))

# ==== 图 1：公开数据中的重疊区 | Fig. 1: the overlap area in the public data ====
ll    <- read_csv(file.path(OUT, "01_ll_long.csv"), show_col_types = FALSE)
share <- read.csv(file.path(OUT, "02_overlap_share_by_year.csv"))
flags <- read.csv(file.path(OUT, "02_flag_ratio_2015_2024.csv"))
cells <- ll |> filter(year >= 2015, year <= 2024) |> group_by(rfmo, lat_sw, lon_sw) |>
  summarise(hooks = sum(hooks), .groups = "drop") |> mutate(lon360 = lon_sw %% 360)
pal <- colorRampPalette(c("#F7FBFF", "#C6DBEF", "#6BAED6", "#2171B5", "#08306B"))(100)
vmax <- max(cells$hooks) / 1e6
colour_of <- round(99 * sqrt(cells$hooks / 1e6 / vmax)) + 1   # 平方根色标 | square-root colour scale
cells$col <- pal[colour_of]

png(file.path(FIG, "Fig1_public_data.png"), width = 7.2, height = 7.6, units = "in", res = 300, pointsize = 7)
layout(matrix(c(1, 1,
                2, 3,
                4, 4), 3, byrow = TRUE), heights = c(1.15, 1, 1.15))
# a) 2015–2024 每个 5° 格的钩数；两组织都有记录的格：左上三角 WCPFC，右下三角 IATTC
# a) hooks per 5° cell, 2015–2024; cells in both products: upper-left triangle WCPFC, lower-right IATTC
par(mar = c(3, 4, 2.5, 6))
plot(NA, xlim = c(140, 290), ylim = c(-55, 5), asp = 1, xaxs = "i", yaxs = "i", axes = FALSE, xlab = "", ylab = "")
wc <- cells[cells$rfmo == "WCPFC", ]; ic <- cells[cells$rfmo == "IATTC", ]
key_w <- paste(wc$lat_sw, wc$lon360); key_i <- paste(ic$lat_sw, ic$lon360)
for (i in seq_len(nrow(wc))) {
  j <- match(key_w[i], key_i)
  if (!is.na(j)) polygon(wc$lon360[i] + c(0, 0, 5), wc$lat_sw[i] + c(0, 5, 5), col = wc$col[i], border = NA)
  else rect(wc$lon360[i], wc$lat_sw[i], wc$lon360[i] + 5, wc$lat_sw[i] + 5, col = wc$col[i], border = NA)
}
for (i in seq_len(nrow(ic))) {
  j <- match(key_i[i], key_w)
  if (!is.na(j)) polygon(ic$lon360[i] + c(0, 5, 5), ic$lat_sw[i] + c(0, 0, 5), col = ic$col[i], border = NA)
  else rect(ic$lon360[i], ic$lat_sw[i], ic$lon360[i] + 5, ic$lat_sw[i] + 5, col = ic$col[i], border = NA)
}
eez_file <- Sys.getenv("EEZ_SHP", "")
if (nzchar(eez_file) && file.exists(eez_file) && requireNamespace("sf", quietly = TRUE)) {
  sf::sf_use_s2(FALSE)                                   # 0–360 经度按平面坐标处理 | treat 0–360 longitudes as planar
  eez <- sf::st_read(eez_file, quiet = TRUE)
  eez <- suppressWarnings(suppressMessages(sf::st_crop(eez, xmin = 140, xmax = 290, ymin = -55, ymax = 5)))
  plot(sf::st_geometry(eez), add = TRUE, border = "#8C8C8C", col = NA, lwd = .3)
}
# 公约区边界：WCPFC 实线、IATTC 虚线，交集为重疊区 | convention boundaries: WCPFC solid, IATTC dashed; intersection = overlap area
lines(c(210, 210, 230, 230), c(0, -4, -4, -55), col = COL[1], lwd = 1.6)
lines(c(210, 210, 290), c(5, -50, -50), col = COL[2], lwd = 1.6, lty = 2)
text(220, -27, "Overlap\narea", font = 2)
text(175, 7.5, "WCPFC Convention Area", col = COL[1], xpd = TRUE); text(250, 7.5, "IATTC Convention Area", col = COL[2], xpd = TRUE)
axis(1, at = c(150, 180, 210, 240, 270), labels = c("150°E", "180°", "150°W", "120°W", "90°W"))
axis(2, at = c(-50, -25, 0), labels = c("50°S", "25°S", "0°"), las = 1); box(bty = "l")
text(grconvertX(.005, "nfc", "user"), grconvertY(.995, "nfc", "user"), "a)", adj = c(0, 1), font = 2, cex = 8 / 7, xpd = NA)
# 色标：地图右侧，平方根刻度 | colour bar to the right of the map, square-root scale
yb <- seq(-55, 5, length.out = 101)
rect(294, yb[-101], 299, yb[-1], col = pal, border = NA, xpd = NA); rect(294, -55, 299, 5, xpd = NA)
ticks <- c(0, 25, 50, 75, 100, 125, 150, 175)
yt <- -55 + 60 * sqrt(ticks / vmax)
segments(299, yt, 300.5, yt, xpd = NA); text(301.5, yt, ticks, adj = 0, xpd = NA)
text(313, -25, "Hooks (millions)", srt = 90, xpd = NA)
# b、c) 重疊区占各自产品的份额（面积加权）| b, c) overlap share of each product (area-weighted)
for (j in 1:2) {
  metric <- c("share_hooks", "share_albmt")[j]
  par(mar = c(4, 4.5, 2.5, 1))
  plot(NA, xlim = c(1990, 2024), ylim = c(0, 100), xlab = "Year", las = 1, xaxs = "i",
       ylab = c("Overlap share of hooks (%)", "Overlap share of albacore weight (%)")[j], bty = "l")
  for (k in 1:2) {
    v <- share[share$rfmo == c("WCPFC", "IATTC")[k] & share$scheme == "area", ]; v <- v[order(v$year), ]
    lines(v$year, 100 * v[[metric]], col = COL[k], lty = k, lwd = 1.6)
  }
  if (j == 1) legend("topright", c("WCPFC", "IATTC"), col = COL[1:2], lty = 1:2, lwd = 1.6, bty = "n")
  text(grconvertX(.005, "nfc", "user"), grconvertY(.995, "nfc", "user"), c("b)", "c)")[j], adj = c(0, 1), font = 2, cex = 8 / 7, xpd = NA)
}
# d) 十年都有正值的五个船旗：IATTC / WCPFC 比值 | d) five flags positive in all ten years: IATTC / WCPFC ratio
order_flags <- c("Chinese Taipei", "French Polynesia", "China", "Korea", "Vanuatu")
par(mar = c(4, 9, 2.5, 1))
plot(NA, xlim = c(.8, 1.9), ylim = c(length(order_flags) + .5, .5), xaxs = "i", yaxs = "i", axes = FALSE,
     xlab = "Ratio of IATTC to WCPFC totals in the overlap area, 2015–2024", ylab = "")
for (i in seq(1, length(order_flags), 2)) rect(.8, i - .5, 1.9, i + .5, col = "#EFEFEF", border = NA)
bw <- .26
metrics_d <- c("hooks", "alb_mt", "alb_n"); col_d <- COL[1:3]
for (m in 1:3) for (i in seq_along(order_flags)) {
  v <- flags$ratio[flags$flag_label == order_flags[i] & flags$metric == metrics_d[m]]
  if (length(v) && !is.na(v)) rect(1, i + (m - 2) * bw - bw / 2, v, i + (m - 2) * bw + bw / 2, col = col_d[m], border = NA)
}
abline(v = 1, lwd = 1)
axis(1, at = c(.8, 1, 1.2, 1.4, 1.6, 1.8)); axis(2, at = seq_along(order_flags), labels = order_flags, las = 1, tick = FALSE)
legend("topright", c("Hooks", "Albacore weight", "Albacore number"), fill = col_d, border = NA, bty = "n")
text(grconvertX(.005, "nfc", "user"), grconvertY(.995, "nfc", "user"), "d)", adj = c(0, 1), font = 2, cex = 8 / 7, xpd = NA)
invisible(dev.off())

# ==== 图 2：模拟设计 | Fig. 2: simulation design ====
png(file.path(FIG, "Fig2_design.png"), width = 7.2, height = 3.9, units = "in", res = 300, pointsize = 7)
par(mar = c(0, 0, 0, 0))
plot(NA, xlim = c(0, 10), ylim = c(0, 7), axes = FALSE, xlab = "", ylab = "")
text(0, 6.85, "a)", font = 2, adj = 0, cex = 8 / 7); text(.4, 6.85, "Three communicating regions and four fleets", adj = 0)
rx <- c(.8, 3.9, 7); fills <- c("#E7F1FA", "#F0EAF7", "#FFF1DA"); edges <- c("#0072B2", "#8064A2", "#A66C00")
titles_r <- c("R1  WCPFC only", "R2  Overlap", "R3  EPO only"); fleets_r <- c("Fleet W", "Fleets OW and OI", "Fleet E")
for (i in 1:3) {
  rect(rx[i], 5.55, rx[i] + 2.25, 6.5, col = fills[i], border = edges[i], lwd = 1.1)
  text(rx[i] + 1.125, 6.17, titles_r[i], font = 2); text(rx[i] + 1.125, 5.85, fleets_r[i])
}
arrows(3.05, 6, 3.9, 6, length = .06, code = 3); arrows(6.15, 6, 7, 6, length = .06, code = 3)
tags <- data.frame(x = c(3.47, 6.57, 4.25, 7.35, 4.55), y = c(6.32, 6.32, 6.66, 6.66, 4.0), s = c("S3", "S3", "S1", "S4", "S2"))
for (i in seq_len(nrow(tags))) {
  rect(tags$x[i] - .17, tags$y[i] - .12, tags$x[i] + .17, tags$y[i] + .12, col = "white", border = "#444444", lwd = .8)
  text(tags$x[i], tags$y[i], tags$s[i], font = 2)
}
text(5, 5.2, "S1  Choice of applicable measures     ·     S3  Time-varying movement     ·     S4  EPO capacity growth")
text(0, 4.7, "b)", font = 2, adj = 0, cex = 8 / 7); text(.4, 4.7, "Annual feedback with management decisions every three years", adj = 0)
boxes <- data.frame(x = c(.6, 6.5, 6.5, .6), y = c(3.25, 3.25, .65, .65),
                    a = c("Operating model", "Observations", "Estimation method", "Management and fishing"),
                    b = c("Survival, recruitment and catch", "Catch and abundance indices", "Age-structured population update",
                          "HCR → TAC → allocation → catch"))
for (i in 1:4) {
  rect(boxes$x[i], boxes$y[i], boxes$x[i] + 2.9, boxes$y[i] + 1, col = "#F2F2F2", border = "#555555", lwd = 1.1)
  text(boxes$x[i] + 1.45, boxes$y[i] + .66, boxes$a[i], font = 2); text(boxes$x[i] + 1.45, boxes$y[i] + .30, boxes$b[i])
}
arrows(c(3.5, 7.95, 6.5, 2.05), c(3.75, 3.25, 1.15, 1.65), c(6.5, 7.95, 3.5, 2.05), c(3.75, 1.65, 1.15, 3.25), length = .08)
text(4.85, 4.0, "Catch misreporting", adj = 0); text(5, 3.52, "Catch and indices")
text(8.1, 2.45, "Reported catch\nand indices", adj = 0); text(5, .9, "Estimated state"); text(2.2, 2.45, "Removals", adj = 0)
invisible(dev.off())

# ==== 图 3：参考情景（C05）第 1 个重复，MP-A 与 MP-B | Fig. 3: replicate 1 of the reference scenario, MP-A and MP-B ====
png(file.path(FIG, "Fig3_baseline.png"), width = 7.2, height = 4.55, units = "in", res = 300, pointsize = 7)
par(mfrow = c(2, 3), mar = c(4, 5, 2.5, 1))
for (i in 1:2) {
  mp <- MPS[i]; d <- traj[traj$case_id == "C05" & traj$mp == mp, ]
  for (j in 1:3) {
    if (j == 1) {               # 耗竭与目标（虚线）| depletion and target (dashed)
      plot(d$year, 100 * d$D, type = "l", col = COL[i], lwd = 2.2, ylim = c(40.2, 46), xlim = c(1, 20), las = 1, bty = "l",
           xlab = "Projection year", ylab = expression(SB / SB[F == 0] ~ "(%)"))
      lines(d$year, 100 * d$target_D, lty = 2, lwd = 1.3)
    } else {                    # 船队渔获（实线）与配额（虚线阶梯）| fleet catch (solid) and quota (dashed steps)
      fl <- c("W", "OW")[j - 1]
      plot(d$year, 100 * d[[paste0("C_", fl)]] / d$last_total_catch, type = "l", col = COL[i], lwd = 2.2, xlim = c(1, 20),
           ylim = c(0, if (fl == "W") 90 else 12), las = 1, bty = "l",
           xlab = "Projection year", ylab = "Catch or quota\n(% of final historical catch)")
      lines(d$year, 100 * d[[paste0("attempted_TAC_", fl)]] / d$last_total_catch, type = "s", lty = 2, lwd = 1.3)
    }
    title(paste(substr(mp, 1, 4), c("Depletion", "W fleet", "OW fleet")[j]), font.main = 1, line = .6)
    text(grconvertX(.005, "nfc", "user"), grconvertY(.995, "nfc", "user"), paste0(letters[(i - 1) * 3 + j], ")"), adj = c(0, 1), font = 2, cex = 8 / 7, xpd = NA)
  }
}
invisible(dev.off())

# ==== 图 4：参考配置下实施因子的边际效应 | Fig. 4: marginal effects of the implementation factors, reference case ====
contrasts <- c("S1+", "S2-", "S2+", "S3+", "S4+")
row_labels <- c("S1 Choice of\napplicable measures", "S2 Under-reporting", "S2 Over-reporting", "S3 Time-varying\nmovement", "S4 EPO capacity growth")
metrics3 <- c("mean_D", "relative_catch", "mae_target")
titles3 <- c("Depletion", "Catch", "Absolute target error")
xlabs3 <- c("Change (percentage points)", "Change (% of final historical catch)", "Change (percentage points)")
d4 <- marg[marg$config_id == "reference", ]
png(file.path(FIG, "Fig4_factor_effects.png"), width = 7.2, height = 3.9, units = "in", res = 300, pointsize = 7)
layout(matrix(1:3, 1), widths = c(1.55, 1, 1)); par(oma = c(0, 0, 2, 0))
for (j in 1:3) {
  par(mar = c(4.5, if (j == 1) 12 else 1, 4, 1))
  v <- d4[d4$metric == metrics3[j], ]
  xr <- range(100 * c(v$estimate - 1.96 * v$mcse, v$estimate + 1.96 * v$mcse, 0))
  plot(NA, xlim = xr + c(-.05, .05) * diff(xr), ylim = c(5.5, .5), yaxs = "i", axes = FALSE, xlab = xlabs3[j], ylab = "")
  for (i in c(2, 4)) rect(par("usr")[1], i - .5, par("usr")[2], i + .5, col = "#EFEFEF", border = NA)
  abline(v = 0, lty = 3, col = "#555555")
  for (k in 1:4) for (i in 1:5) {
    r <- v[v$mp == MPS[k] & v$contrast == contrasts[i], ]
    y <- i + (k - 2.5) * .17
    segments(100 * (r$estimate - 1.96 * r$mcse), y, 100 * (r$estimate + 1.96 * r$mcse), y, col = COL[k])
    points(100 * r$estimate, y, pch = PCH[k], col = COL[k], cex = 1.3)
  }
  axis(1); box(bty = "l")
  if (j == 1) axis(2, at = 1:5, labels = row_labels, las = 1, tick = FALSE)
  title(titles3[j], font.main = 1, line = .6)
  text(grconvertX(.005, "nfc", "user"), grconvertY(.995, "nfc", "user"), paste0(letters[j], ")"), adj = c(0, 1), font = 2, cex = 8 / 7, xpd = NA)
}
legend(grconvertX(.5, "ndc", "user"), grconvertY(1, "ndc", "user"), LAB, col = COL, pch = PCH, horiz = TRUE,
       bty = "n", xjust = .5, yjust = 1, xpd = NA)                 # 图例在上方外边距 | legend in the top outer margin
invisible(dev.off())

# ==== 图 5：MP-B 只开误报 | Fig. 5: MP-B with misreporting only ====
d5 <- absd[absd$config_id == "reference" & absd$mp == MPS[2] & absd$S1 == 0 & absd$S3 == 0 & absd$S4 == 0, ]
xl <- c("Under-reporting", "Unbiased\nreporting", "Over-reporting")
png(file.path(FIG, "Fig5_reporting.png"), width = 7.2, height = 4.9, units = "in", res = 300, pointsize = 7)
par(mfrow = c(2, 2), mar = c(4.5, 5, 2.5, 1))
panels5 <- data.frame(metric = c("mean_D", "relative_catch", "mae_target"), title = c("Depletion", "Catch", "Absolute target error"),
                      ylab = c("SB/SB(F=0) (%)", "% of final historical catch", "Percentage points"),
                      lo = c(41, 99, 0), hi = c(45.6, 108.5, 1.9), digits = c(1, 1, 2))
for (i in 1:3) {
  v <- d5[d5$metric == panels5$metric[i], ]; v <- v[order(v$S2), ]
  y <- 100 * v$estimate; e <- 100 * 1.96 * v$mcse
  plot(1:3, y, type = "l", col = "#9A9A9A", lwd = 2, xlim = c(.55, 3.45), ylim = c(panels5$lo[i], panels5$hi[i]),
       axes = FALSE, xlab = if (i == 3) "Reported overlap catch" else "", ylab = panels5$ylab[i])
  arrows(1:3, y - e, 1:3, y + e, angle = 90, code = 3, length = .02, col = COL[2])
  points(1:3, y, pch = 15, col = COL[2], cex = 1.5)
  text(1:3, y, paste0(formatC(y, format = "f", digits = panels5$digits[i]), if (i < 3) "%" else ""), pos = 3, offset = .7)
  if (i == 1) { abline(h = 43.2, lty = 2, lwd = 1.4); text(3.4, 43.4, "Target 43.2%", adj = 1) }
  if (i == 2) { abline(h = 100, lty = 2, lwd = 1.4); text(3.4, 100.4, "Historical catch", adj = 1) }
  axis(1, at = 1:3, labels = xl, padj = .5); axis(2, las = 1); box(bty = "l")
  title(panels5$title[i], font.main = 1, line = .6)
  text(grconvertX(.005, "nfc", "user"), grconvertY(.995, "nfc", "user"), paste0(letters[i], ")"), adj = c(0, 1), font = 2, cex = 8 / 7, xpd = NA)
}
qW <- d5[d5$metric == "quota_fully_caught_W", ]; qW <- qW[order(qW$S2), ]
qO <- d5[d5$metric == "quota_fully_caught_OW", ]; qO <- qO[order(qO$S2), ]
bars <- rbind(100 * qW$estimate, 100 * qO$estimate)
bp <- barplot(bars, beside = TRUE, col = c("#4D4D4D", "#B8B8B8"), border = NA, ylim = c(0, 128), axes = FALSE,
              axisnames = FALSE, xlab = "Reported overlap catch", ylab = "% of evaluation years")
axis(1, at = colMeans(bp), labels = xl, padj = .5, tick = FALSE)
text(bp, bars, paste0(round(bars), "%"), pos = 3)
axis(2, at = seq(0, 100, 20), las = 1)
legend("top", c("W fleet", "OW fleet"), fill = c("#4D4D4D", "#B8B8B8"), border = NA, horiz = TRUE, bty = "n")
title("Years with quota fully caught", font.main = 1, line = .6)
text(grconvertX(.005, "nfc", "user"), grconvertY(.995, "nfc", "user"), "d)", adj = c(0, 1), font = 2, cex = 8 / 7, xpd = NA)
invisible(dev.off())

# ==== 图 6：合并情景（C23）在 12 个配置下的效应 | Fig. 6: the combined scenario (C23) across the 12 configurations ====
cfg_rows <- CONFIGS$config_id
cfg_labels <- c("Reference case", "Adult M = 0.24 per year", "Adult M = 0.36 per year", "Adult M = 0.54 per year",
                "Steepness h = 0.65", "Steepness h = 0.95", "Connectivity × 0.5", "Connectivity × 2",
                "Selectivity a50 = 2.5", "Selectivity a50 = 4.5", "Errors in history\nand projection", "Movement information\nin estimation")
d6 <- paired[paired$case_id == "C23", ]
png(file.path(FIG, "Fig6_sensitivity.png"), width = 7.2, height = 6.25, units = "in", res = 300, pointsize = 7)
layout(matrix(1:2, 1), widths = c(1.6, 1)); par(oma = c(0, 0, 2, 0))
for (j in 1:2) {
  par(mar = c(4.5, if (j == 1) 12 else 1, 4, 1))
  v <- d6[d6$metric == metrics3[j], ]
  xr <- range(100 * c(v$estimate - 1.96 * v$mcse, v$estimate + 1.96 * v$mcse, 0))
  plot(NA, xlim = xr + c(-.05, .05) * diff(xr), ylim = c(12.5, .5), yaxs = "i", axes = FALSE, xlab = xlabs3[j], ylab = "")
  for (i in seq(2, 12, 2)) rect(par("usr")[1], i - .5, par("usr")[2], i + .5, col = "#EFEFEF", border = NA)
  abline(v = 0, lty = 3, col = "#555555"); abline(h = 6.5, col = "#999999")   # 生物学配置 / 其他敏感性 | biology / other sensitivities
  for (k in 1:4) for (i in 1:12) {
    r <- v[v$mp == MPS[k] & v$config_id == cfg_rows[i], ]
    y <- i + (k - 2.5) * .16
    segments(100 * (r$estimate - 1.96 * r$mcse), y, 100 * (r$estimate + 1.96 * r$mcse), y, col = COL[k])
    points(100 * r$estimate, y, pch = PCH[k], col = COL[k], cex = 1.2)
  }
  axis(1); box(bty = "l")
  if (j == 1) axis(2, at = 1:12, labels = cfg_labels, las = 1, tick = FALSE)
  title(titles3[j], font.main = 1, line = .6)
  text(grconvertX(.005, "nfc", "user"), grconvertY(.995, "nfc", "user"), paste0(letters[j], ")"), adj = c(0, 1), font = 2, cex = 8 / 7, xpd = NA)
}
legend(grconvertX(.5, "ndc", "user"), grconvertY(1, "ndc", "user"), LAB, col = COL, pch = PCH, horiz = TRUE,
       bty = "n", xjust = .5, yjust = 1, xpd = NA)                 # 图例在上方外边距 | legend in the top outer margin
invisible(dev.off())
cat("12 done:", FIG, "\n")
