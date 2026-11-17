# 09_fixed_shares.R — 补充实验：MP-B 固定上报份额（补充说明 6、补充表 13）
# 09_fixed_shares.R — supplementary experiment: MP-B with fixed reported shares (Supplementary Note 6, Table 13)
#
# MP-B 每 3 年按最近 3 年上报渔获份额重新分配共同 TAC。这里把份额固定在预测开始时的值（历史最后 3 年），
# 其余（HCR、TAC 换算、重疊区按容量分配、误报进入估计方法）不变。需要先运行 07（CONFIG=reference）。
# MP-B re-divides the joint TAC every 3 years by the last 3 years' reported shares. Here the shares are frozen at
# their values at the start of the projection (last 3 historical years); everything else is unchanged.
# Needs step 07 for CONFIG=reference first.
#   N_REP=300 Rscript R/09_fixed_shares.R      # 约 5 分钟 | about 5 min

if (!exists("N_REP")) N_REP <- as.integer(Sys.getenv("N_REP", "300"))
CONFIG <- "reference"
RUN_CELLS <- NULL
source("R/00_settings.R")
RUN_CELLS <- design$case_id                       # 全部 24 个情景 | all 24 scenarios
RUN_MPS <- "MP-B joint TAC"                       # 只跑 MP-B | MP-B only
FIXED_SHARES <- TRUE                              # 06 中 MP-B 的分配改用固定份额 | step 06 uses the fixed shares for MP-B

# ---- 1. 运行 03–06 | run steps 03–06 ----
source("R/03_parameters.R"); source("R/04_initialise.R"); source("R/05_history.R")
cat("fixed shares of replicate 1 (W, OW, OI, E):",
    round(colSums(init[[1]]$catch_history) / sum(init[[1]]$catch_history), 4), "\n")
source("R/06_closed_loop.R")
fixed <- metrics
dest <- file.path(OUT, "09_fixed_shares"); dir.create(dest, showWarnings = FALSE)
con <- gzfile(file.path(dest, "replicate_metrics.csv.gz"), "wt"); write.csv(fixed, con, row.names = FALSE, na = ""); close(con)

# ---- 2. 滚动份额 vs 固定份额的边际效应（补充表 13）| marginal effects, rolling vs fixed shares (Supplementary Table 13) ----
rolling <- read.csv(gzfile(file.path(OUT, "07_reference", "replicate_metrics.csv.gz"))) |>
  filter(mp == "MP-B joint TAC", rep <= N_REP)
stopifnot(max(rolling$rep) == N_REP)
effects <- list()
for (variant in c("rolling", "fixed")) {
  m <- if (variant == "rolling") rolling else fixed
  for (f in c("S1", "S2", "S3", "S4")) for (o in c("mean_D", "relative_catch", "mae_target")) {
    keys <- c("rep", setdiff(c("S1", "S2", "S3", "S4"), f))
    base <- m[m[[f]] == 0, c(keys, o)]; names(base)[ncol(base)] <- "base"
    z <- inner_join(m[m[[f]] != 0, c(keys, f, o)], base, by = keys)
    z$delta <- z[[o]] - z$base
    effects[[length(effects) + 1]] <- z |> group_by(rep, level = .data[[f]]) |> summarise(delta = mean(delta), .groups = "drop") |>
      group_by(level) |> summarise(estimate = mean(delta), mcse = sd(delta) / sqrt(n()), .groups = "drop") |>
      mutate(variant = variant, factor = paste0(f, ifelse(level > 0, "+", "-")), metric = o)
  }
}
# 百分点或 % 历史末年渔获；小于 1e-12 的数值误差记为 0 | percentage points or % of final historical catch; |x| < 1e-12 shown as 0
effects <- bind_rows(effects) |> mutate(value = ifelse(abs(estimate) < 1e-12, 0, signif(100 * estimate, 3)))
write.csv(effects, file.path(dest, "effects_rolling_vs_fixed.csv"), row.names = FALSE)
print(as.data.frame(effects |> select(factor, metric, variant, value) |> pivot_wider(names_from = variant, values_from = value)))
cat("09 done:", dest, "\n")
