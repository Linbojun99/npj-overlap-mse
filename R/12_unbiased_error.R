# 12_unbiased_error.R — 补充实验：重疊区上报渔获的无偏误差（补充说明 5、补充表 11、图 4、图 5）
# 12_unbiased_error.R — supplementary experiment: unbiased error in the reported overlap catch
#                       (Supplementary Note 5, Supplementary Table 11, Figs 4 and 5)
#
# 主实验的少报与多报是单向误差：b = m (0.5 + u)，m = −0.25 或 +0.25，u 为 0–1 的均匀随机数。这里加一个无偏水平：
# 同一个 u 换成 3u − 2，b = 0.25 (0.5 + 3u − 2) = 0.75u − 0.375，在 −0.375 与 0.375 之间均匀分布、均值为零，
# 正好是两个单向水平合起来的范围。只有上报渔获改变，其余随机数不变，所以每个重复都与主实验中准确上报的同一重复配对。
# 参考配置；S1、S3、S4 的 8 种组合 × 4 个候选 MP × N_REP 个重复（论文 300，共 9,600 条轨迹）。
# 需要先运行 07（CONFIG=reference，N_REP 相同）。
# Under- and over-reporting in the main experiment are errors in one direction: b = m (0.5 + u), m = −0.25 or
# +0.25, u uniform on 0–1. This adds an unbiased level: the same u is replaced by 3u − 2, so
# b = 0.25 (0.5 + 3u − 2) = 0.75u − 0.375, uniform between −0.375 and 0.375 with mean zero, the range spanned by
# the two directional levels. Only the reported catch changes and every other random number is the same, so each
# replicate is paired with the same replicate of the main experiment with accurate reporting.
# Reference case; the 8 combinations of S1, S3 and S4 × 4 candidate MPs × N_REP replicates (300 in the Article,
# 9,600 trajectories). Needs step 07 for CONFIG=reference (same N_REP) first.
#   N_REP=300 Rscript R/12_unbiased_error.R      # 约 2 分钟 | about 2 min

if (!exists("N_REP")) N_REP <- as.integer(Sys.getenv("N_REP", "300"))
CONFIG <- "reference"
RUN_CELLS <- NULL
source("R/00_settings.R")
RUN_CELLS <- design$case_id[design$S2 == 0]           # 准确上报的 8 个情景 C05–C08、C17–C20 | the 8 accurate-reporting scenarios
RUN_MPS <- c("MP-A WCPFC TAC", "MP-B joint TAC", "MP-C TAC with deduction", "MP-D overlap effort cap")
FIXED_SHARES <- FALSE                                 # MP-B 用按上报渔获更新的份额 | MP-B uses shares updated from reported catch
UNBIASED_ERROR <- TRUE                                # 06 把这 8 个情景的上报渔获乘 (1 + b) | step 06 multiplies their reported overlap catch by (1 + b)

# ---- 1. 运行 03–05，检查误差的范围 | run steps 03–05 and check the range of the error ----
source("R/03_parameters.R"); source("R/04_initialise.R"); source("R/05_history.R")
proj <- n_hist + seq_len(n_proj)                      # 预测年在随机数中的位置 | positions of the projection years in the draws
b <- matrix(NA_real_, N_REP, n_proj)                  # 每个重复、每个预测年的 b（与 06 的算式相同）| b by replicate and year (as in step 06)
b_mean <- numeric(N_REP)                              # 每个重复第 1–20 年 b 的均值 | mean b of each replicate over years 1–20
for (r in seq_len(N_REP)) {
  b[r, ] <- abs(misreporting) * (.5 + (3 * draws[[r]]$reporting_u[proj] - 2))
  b_mean[r] <- mean(b[r, ])
}
stopifnot(all(abs(b) <= .375 + 1e-12))                # 每年都在 ±0.375 之内 | inside ±0.375 in every year
# b 在所有年与重复上的均值，以及第 1 个重复的范围 | mean b over years and replicates, and the range in replicate 1
summary_b <- data.frame(mean_b = mean(b_mean), replicate1_min = min(b[1, ]), replicate1_max = max(b[1, ]))
print(summary_b)

# ---- 2. 闭环（06）| closed loop (step 06) ----
source("R/06_closed_loop.R")
unbiased <- metrics
unbiased$S2 <- "unbiased_error"                       # S2 的第四个水平 | the fourth level of S2
unbiased$mean_reporting_error <- b_mean[unbiased$rep]  # 该重复第 1–20 年 b 的均值 | mean b of the replicate over years 1–20
dest <- file.path(OUT, "12_unbiased_error"); dir.create(dest, showWarnings = FALSE)
con <- gzfile(file.path(dest, "replicate_metrics.csv.gz"), "wt"); write.csv(unbiased, con, row.names = FALSE, na = ""); close(con)
write.csv(trajectories, file.path(dest, "first_replicate_trajectories.csv"), row.names = FALSE, na = "")
write.csv(diagnostics, file.path(dest, "diagnostics.csv"), row.names = FALSE, na = "")
write.csv(summary_b, file.path(dest, "reporting_error_summary.csv"), row.names = FALSE)

# ---- 3. 与主实验配对：三个误差水平相对准确上报的效应（图 4、补充表 11）----
# ---- 3. pairing with the main experiment: effects of the three error levels against accurate reporting (Fig. 4, Supplementary Table 11) ----
# 同一 MP、重复与 S1、S3、S4 组合内相减，先在重复内对 8 种组合平均，再跨重复求均值与 MCSE；
# 少报与多报的效应就是 08 中参考配置的 S2 效应。
# Differences within MP, replicate and combination of S1, S3 and S4 are averaged over the 8 combinations within
# replicate, then the mean and MCSE are taken across replicates; the effects of under- and over-reporting are
# the S2 effects of step 08 in the reference case.
main <- read.csv(gzfile(file.path(OUT, "07_reference", "replicate_metrics.csv.gz"))) |> filter(rep <= N_REP)
stopifnot(max(main$rep) == N_REP)
outcomes <- c("mean_D", "relative_catch", "quota_fully_caught_W", "quota_fully_caught_OW", "uncaught_quota_ratio")
keys <- c("mp", "rep", "S1", "S3", "S4")
accurate <- main[main$S2 == 0, ]
error_levels <- list(unbiased_error = unbiased, under_reporting = main[main$S2 == -1, ], over_reporting = main[main$S2 == 1, ])
effects <- list()
for (lev in names(error_levels)) for (o in outcomes) {
  base <- accurate[c(keys, o)]; names(base)[ncol(base)] <- "base"
  z <- inner_join(error_levels[[lev]][c(keys, o)], base, by = keys)
  stopifnot(nrow(z) == 8 * length(RUN_MPS) * N_REP)  # 每个 MP、重复各 8 对 | 8 pairs per MP and replicate
  z$delta <- z[[o]] - z$base
  effects[[paste(lev, o)]] <- z |>
    group_by(mp, rep) |> summarise(delta = mean(delta), .groups = "drop") |>       # 重复内 | within replicate
    group_by(mp) |> summarise(estimate = mean(delta), n_rep = n(), mcse = sd(delta) / sqrt(n()), .groups = "drop") |>
    mutate(level = lev, metric = o)
}
effects <- bind_rows(effects)
write.csv(effects, file.path(dest, "effects.csv"), row.names = FALSE)

# ---- 4. S1、S3、S4 都关时四个上报水平的均值（图 5、补充表 11）----
# ---- 4. means of the four reporting levels with S1, S3 and S4 off (Fig. 5, Supplementary Table 11) ----
off <- list()
for (lev in c("accurate", names(error_levels))) {
  m <- if (lev == "accurate") accurate else error_levels[[lev]]
  off[[lev]] <- m[m$S1 == 0 & m$S3 == 0 & m$S4 == 0, c("mp", outcomes)] |> mutate(level = lev)
}
levels_off <- bind_rows(off) |> pivot_longer(all_of(outcomes), names_to = "metric") |>
  group_by(level, mp, metric) |> summarise(mean = mean(value), n_rep = n(), mcse = sd(value) / sqrt(n()), .groups = "drop")
write.csv(levels_off, file.path(dest, "levels_other_factors_off.csv"), row.names = FALSE)

# ---- 5. 补充表 11 | Supplementary Table 11 ----
# 前两行是效应（均值 ± MCSE；百分点或 % 历史末年渔获），其余是 S1、S3、S4 都关时第 11–20 年的均值（%）
# The first two rows are effects (mean ± MCSE; percentage points or % of final historical catch), the others are
# means over years 11–20 with S1, S3 and S4 off (%).
rows11 <- c(mean_D = "Change in SB/SB_F=0 (percentage points)", relative_catch = "Change in catch (% of final historical catch)")
t11a <- effects |> filter(level == "unbiased_error", metric %in% names(rows11)) |>
  mutate(value = sprintf("%+.3f", 100 * estimate), value = ifelse(as.numeric(value) == 0, "0.000", value),   # 0 不带正负号 | zero unsigned
         se = sprintf("%.3f", 100 * mcse), se = ifelse(as.numeric(se) == 0, "<0.001", se),
         outcome = rows11[metric], cell = paste(value, "±", se))
rows11b <- c(mean_D = "SB/SB_F=0 (%)", relative_catch = "Catch (% of final historical catch)",
             quota_fully_caught_W = "Years in which fleet W takes its full quota (%)",
             quota_fully_caught_OW = "Years in which fleet OW takes its full quota (%)",
             uncaught_quota_ratio = "Uncaught part of the TAC (%)")
t11b <- levels_off |> filter(level == "unbiased_error") |>
  mutate(outcome = rows11b[metric], cell = ifelse(metric == "uncaught_quota_ratio", sprintf("%.2f", 100 * mean), sprintf("%.1f", 100 * mean)))
table11 <- bind_rows(t11a, t11b) |> transmute(outcome = factor(outcome, c(rows11, rows11b)), mp = substr(mp, 1, 4), cell) |>
  pivot_wider(names_from = mp, values_from = cell) |> arrange(outcome)
write.csv(table11, file.path(dest, "supplementary_table_11.csv"), row.names = FALSE)
cat("\n== Supplementary Table 11: unbiased reporting error ==\n")
print(as.data.frame(table11), right = FALSE)

# ---- 6. 图 5 的数值：MP-B，S1、S3、S4 都关 | values of Fig. 5: MP-B with S1, S3 and S4 off ----
cat("\n== MP-B with S1, S3 and S4 off (percent) ==\n")
print(as.data.frame(levels_off |> filter(mp == "MP-B joint TAC") |> mutate(value = round(100 * mean, 2)) |>
  select(level, metric, value) |> pivot_wider(names_from = metric, values_from = value)))
cat("12 done:", dest, "\n")
