# 13_figures_supplementary.R — 补充图 1–5 | Supplementary Figs 1–5
# 需要 07（reference）、08、10（六种渔具结构）与 11 的结果；缺少的图会跳过
# needs the results of steps 07 (reference), 08, 10 (six gear structures) and 11; missing inputs skip that figure

source("R/00_settings.R")
COL  <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7")      # MP-A 蓝、MP-B 橙、MP-C 绿、MP-D 粉 | blue, orange, green, pink
PCH  <- c(16, 15, 17, 18)
MPS  <- c("MP-A WCPFC TAC", "MP-B joint TAC", "MP-C TAC with deduction", "MP-D overlap effort cap")
LAB  <- c("MP-A", "MP-B", "MP-C", "MP-D")
traj   <- read.csv(file.path(OUT, "07_reference", "first_replicate_trajectories.csv"))
paired <- read.csv(file.path(OUT, "08_analysis", "paired.csv"))
info   <- read.csv(file.path(OUT, "08_analysis", "information.csv"))

# ==== 补充图 1：参考情景第 1 个重复，MP-C 与 MP-D | Supp. Fig. 1: replicate 1 of the reference scenario, MP-C and MP-D ====
png(file.path(FIG, "FigS1_baseline_CD.png"), width = 7.2, height = 4.55, units = "in", res = 300, pointsize = 7)
par(mfrow = c(2, 3), mar = c(4, 5, 2.5, 1))
for (i in 3:4) {
  mp <- MPS[i]; d <- traj[traj$case_id == "C05" & traj$mp == mp, ]
  for (j in 1:3) {
    if (j == 1) {
      plot(d$year, 100 * d$D, type = "l", col = COL[i], lwd = 2.2, ylim = c(40.2, 46), xlim = c(1, 20), las = 1, bty = "l",
           xlab = "Projection year", ylab = expression(SB / SB[F == 0] ~ "(%)"))
      lines(d$year, 100 * d$target_D, lty = 2, lwd = 1.3)
    } else {
      fl <- c("W", "OW")[j - 1]
      plot(d$year, 100 * d[[paste0("C_", fl)]] / d$last_total_catch, type = "l", col = COL[i], lwd = 2.2, xlim = c(1, 20),
           ylim = c(0, if (fl == "W") 90 else 12), las = 1, bty = "l",
           xlab = "Projection year", ylab = "Catch or quota\n(% of final historical catch)")
      lines(d$year, 100 * d[[paste0("attempted_TAC_", fl)]] / d$last_total_catch, type = "s", lty = 2, lwd = 1.3)
    }
    title(paste(substr(mp, 1, 4), c("Depletion", "W fleet", "OW fleet")[j]), font.main = 1, line = .6)
    text(grconvertX(.005, "nfc", "user"), grconvertY(.995, "nfc", "user"), paste0(letters[(i - 3) * 3 + j], ")"), adj = c(0, 1), font = 2, cex = 8 / 7, xpd = NA)
  }
}
invisible(dev.off())

# ==== 补充图 2：参考配置下 24 个情景相对 C05 的配对效应 | Supp. Fig. 2: paired effects of the 24 scenarios against C05 ====
d2 <- paired[paired$config_id == "reference", ]
cases <- unique(d2[c("case_id", "S1", "S2", "S3", "S4")])
cases <- cases[order(cases$S3, cases$S4, cases$S2, cases$S1), ]
rep_word <- c("-1" = "under-reporting", "0" = "unbiased", "1" = "over-reporting")
case_labels <- paste0("S1 ", ifelse(cases$S1 == 1, "on", "off"), ", S2 ", rep_word[as.character(cases$S2)],
                      ", S3 ", ifelse(cases$S3 == 1, "on", "off"), ", S4 ", ifelse(cases$S4 == 1, "on", "off"),
                      ifelse(cases$case_id == "C05", " (reference scenario)", ""))
div <- colorRampPalette(c("#B2182B", "#EF8A62", "#FDDBC7", "#F7F7F7", "#D1E5F0", "#67A9CF", "#2166AC"))(101)   # 红–白–蓝 | red–white–blue
png(file.path(FIG, "FigS2_full_factorial.png"), width = 7.2, height = 6.4, units = "in", res = 300, pointsize = 7)
layout(matrix(c(1, 2, 3, 4), 2, byrow = TRUE), widths = c(2.05, 1), heights = c(1, .09))
for (j in 1:2) {
  metric <- c("mean_D", "relative_catch")[j]
  mat_v <- matrix(NA, 24, 4)
  for (i in 1:24) for (k in 1:4)
    mat_v[i, k] <- 100 * d2$estimate[d2$case_id == cases$case_id[i] & d2$mp == MPS[k] & d2$metric == metric]
  lim <- max(abs(mat_v))
  par(mar = c(2.5, if (j == 1) 24 else 1, 3, 1))
  image(1:4, 1:24, t(mat_v[24:1, ]), col = div, zlim = c(-lim, lim), axes = FALSE, xlab = "", ylab = "")
  axis(1, at = 1:4, labels = LAB, tick = FALSE)
  if (j == 1) axis(2, at = 24:1, labels = case_labels, las = 1, tick = FALSE)
  box(); title(c("Depletion", "Catch")[j], font.main = 1, line = .6)
  text(grconvertX(.005, "nfc", "user"), grconvertY(.995, "nfc", "user"), paste0(letters[j], ")"), adj = c(0, 1), font = 2, cex = 8 / 7, xpd = NA)
  assign(paste0("lim", j), lim)
}
for (j in 1:2) {                                         # 横向色标 | horizontal colour bars
  lim <- get(paste0("lim", j))
  par(mar = c(3, if (j == 1) 24 else 1, .5, 1))
  image(seq(-lim, lim, length.out = 101), 1, matrix(1:101), col = div, axes = FALSE, xlab = "", ylab = "")
  axis(1); box()
  mtext(c("Change (percentage points)", "Change (% of final historical catch)")[j], side = 1, line = 2)
}
invisible(dev.off())

# ==== 补充图 3：估计方法得到真实移动信息（S3 开）| Supp. Fig. 3: realized movement supplied to the estimation (S3 on) ====
png(file.path(FIG, "FigS3_movement_information.png"), width = 7.2, height = 2.7, units = "in", res = 300, pointsize = 7)
par(mfrow = c(1, 2), mar = c(4.5, 5, 2.5, 1))
for (j in 1:2) {
  metric <- c("mean_D", "estimation_rmse_fixed")[j]
  v <- info[info$S3 == 1 & info$metric == metric, ]; v <- v[match(MPS, v$mp), ]
  xr <- range(100 * c(v$estimate - 1.96 * v$mcse, v$estimate + 1.96 * v$mcse, 0))
  plot(NA, xlim = xr + c(-.08, .08) * diff(xr), ylim = c(4.6, .4), axes = FALSE, ylab = "",
       xlab = "Observed minus mean movement (percentage points)")
  abline(v = 0, lty = 3)
  segments(100 * (v$estimate - 1.96 * v$mcse), 1:4, 100 * (v$estimate + 1.96 * v$mcse), 1:4, col = COL)
  points(100 * v$estimate, 1:4, pch = PCH, col = COL, cex = 1.2)
  axis(1); axis(2, at = 1:4, labels = LAB, las = 1, tick = FALSE); box(bty = "l")
  title(c("Depletion", "Estimation error (RMSE)")[j], font.main = 1, line = .6)
  text(grconvertX(.005, "nfc", "user"), grconvertY(.995, "nfc", "user"), paste0(letters[j], ")"), adj = c(0, 1), font = 2, cex = 8 / 7, xpd = NA)
}
invisible(dev.off())

# ==== 补充图 4：分开延绳钓与曳绳钓时合并情景（C23 vs C05）的效应 | Supp. Fig. 4: combined scenario (C23 vs C05) with longline and troll ====
structures <- c("reference", "shared_2pct", "shared_5pct", "troll11_2pct", "troll11_5pct", "troll12_2pct", "troll12_5pct")
struct_labels <- c("Four fleets without\na troll component", "Shared selectivity\n2% troll share", "Shared selectivity\n5% troll share",
                   "Troll curve 11\n2% troll share", "Troll curve 11\n5% troll share", "Troll curve 12\n2% troll share",
                   "Troll curve 12\n5% troll share")
files4 <- c(file.path(OUT, "07_reference", "replicate_metrics.csv.gz"),
            file.path(OUT, paste0("10_", structures[-1]), "replicate_metrics.csv.gz"))
if (all(file.exists(files4))) {
  eff4 <- list()
  for (s in seq_along(structures)) {
    m <- read.csv(gzfile(files4[s]))
    a <- m[m$case_id == "C05", c("mp", "rep", "mean_D", "relative_catch", "mae_target")]
    b <- m[m$case_id == "C23", c("mp", "rep", "mean_D", "relative_catch", "mae_target")]
    z <- inner_join(b, a, by = c("mp", "rep"), suffix = c("", "_ref"))
    for (o in c("mean_D", "relative_catch", "mae_target")) {
      z$delta <- z[[o]] - z[[paste0(o, "_ref")]]
      eff4[[paste(s, o)]] <- z |> group_by(mp) |> summarise(estimate = mean(delta), mcse = sd(delta) / sqrt(n()), .groups = "drop") |>
        mutate(structure = structures[s], metric = o)
    }
  }
  eff4 <- bind_rows(eff4)
  write.csv(eff4, file.path(OUT, "13_gear_combined_effects.csv"), row.names = FALSE)
  png(file.path(FIG, "FigS4_gear_structure.png"), width = 7.2, height = 4.4, units = "in", res = 300, pointsize = 7)
  layout(matrix(1:3, 1), widths = c(1.55, 1, 1)); par(oma = c(0, 0, 2, 0))
  for (j in 1:3) {
    metric <- c("mean_D", "relative_catch", "mae_target")[j]
    v <- eff4[eff4$metric == metric, ]
    xr <- range(100 * c(v$estimate - 1.96 * v$mcse, v$estimate + 1.96 * v$mcse, 0))
    par(mar = c(4.5, if (j == 1) 10 else 1, 4, 1))
    plot(NA, xlim = xr + c(-.05, .05) * diff(xr), ylim = c(7.5, .5), yaxs = "i", axes = FALSE, ylab = "",
         xlab = c("Change (percentage points)", "Change (% of final historical catch)", "Change (percentage points)")[j])
    for (i in seq(2, 7, 2)) rect(par("usr")[1], i - .5, par("usr")[2], i + .5, col = "#EFEFEF", border = NA)
    abline(v = 0, lty = 3, col = "#555555")
    for (k in 1:4) for (i in 1:7) {
      r <- v[v$mp == MPS[k] & v$structure == structures[i], ]
      y <- i + (k - 2.5) * .17
      segments(100 * (r$estimate - 1.96 * r$mcse), y, 100 * (r$estimate + 1.96 * r$mcse), y, col = COL[k])
      points(100 * r$estimate, y, pch = PCH[k], col = COL[k], cex = 1.3)
    }
    axis(1); box(bty = "l")
    if (j == 1) axis(2, at = 1:7, labels = struct_labels, las = 1, tick = FALSE)
    title(c("Depletion", "Catch", "Absolute target error")[j], font.main = 1, line = .6)
    text(grconvertX(.005, "nfc", "user"), grconvertY(.995, "nfc", "user"), paste0(letters[j], ")"), adj = c(0, 1), font = 2, cex = 8 / 7, xpd = NA)
  }
  legend(grconvertX(.5, "ndc", "user"), grconvertY(1, "ndc", "user"), LAB, col = COL, pch = PCH, horiz = TRUE,
         bty = "n", xjust = .5, yjust = 1, xpd = NA)
  invisible(dev.off())
} else cat("Supplementary Fig. 4 skipped: run step 10 for all six structures\n")

# ==== 补充图 5：MP 范围外渔获，两侧对称 | Supp. Fig. 5: catch outside the MP's scope, both sides ====
f5 <- file.path(OUT, "11_external_catch", "marginal_effects.csv")
if (file.exists(f5)) {
  ed <- read.csv(f5)
  levels5 <- data.frame(E = c(-1, 1, 0, 0, -1, 1), T = c(0, 0, -1, 1, -1, 1),
                        label = c("EPO outside the overlap area\n(fleet E), −3% per year", "EPO outside the overlap area\n(fleet E), +3% per year",
                                  "WCPFC fishery outside the MP\n(fleet T), −3% per year", "WCPFC fishery outside the MP\n(fleet T), +3% per year",
                                  "Both fleets, −3% per year", "Both fleets, +3% per year"))
  png(file.path(FIG, "FigS5_external_catch.png"), width = 7.2, height = 3.9, units = "in", res = 300, pointsize = 7)
  layout(matrix(1:3, 1), widths = c(1.6, 1, 1)); par(oma = c(0, 0, 2, 0))
  for (j in 1:3) {
    metric <- c("mean_D", "relative_catch_MP_fleets", "relative_catch")[j]
    v <- ed[ed$metric == metric, ]
    xr <- range(100 * c(v$effect - 1.96 * v$mcse, v$effect + 1.96 * v$mcse, 0))
    par(mar = c(4.5, if (j == 1) 12 else 1, 3, 1))
    plot(NA, xlim = xr + c(-.05, .05) * diff(xr), ylim = c(6.5, .5), yaxs = "i", axes = FALSE, ylab = "",
         xlab = c("Change (percentage points)", "Change (% of final historical catch)", "Change (% of final historical catch)")[j])
    for (i in c(1, 3, 5)) rect(par("usr")[1], i - .5, par("usr")[2], i + .5, col = "#F0F0F0", border = NA)
    abline(v = 0)
    for (k in 1:4) for (i in 1:6) {
      r <- v[v$mp == LAB[k] & v$E_level == levels5$E[i] & v$T_level == levels5$T[i], ]
      y <- i + (k - 2.5) * .17
      segments(100 * (r$effect - 1.96 * r$mcse), y, 100 * (r$effect + 1.96 * r$mcse), y, col = COL[k])
      points(100 * r$effect, y, pch = PCH[k], col = COL[k], cex = 1.3)
    }
    axis(1); box(bty = "l")
    if (j == 1) axis(2, at = 1:6, labels = levels5$label, las = 1, tick = FALSE)
    title(c("Depletion", "Catch of the MP fleets", "Total catch")[j], font.main = 1, line = .6)
    text(grconvertX(.005, "nfc", "user"), grconvertY(.995, "nfc", "user"), paste0(letters[j], ")"), adj = c(0, 1), font = 2, cex = 8 / 7, xpd = NA)
  }
  invisible(dev.off())
} else cat("Supplementary Fig. 5 skipped: run step 11 first\n")
cat("13 done:", FIG, "\n")
