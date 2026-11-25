# 08_effects.R — 配对效应、边际效应、补充表 9 与 14，以及正文中的数字（图 4、5，补充图 2–4，补充表 7–10、14、15）
# 08_effects.R — paired effects, marginal effects, Supplementary Tables 9 and 14, and the numbers quoted in the Article
#                (Figs 4 and 5, Supplementary Figs 2–4, Supplementary Tables 7–10, 14 and 15)
#
# 需要 07 对 12 个配置各运行一次（N_REP 相同）| needs step 07 run once for each of the 12 configurations (same N_REP)
# 配对：同一 MP、同一重复内比较（随机数相同）；效应先在重复内平均，再跨重复求均值与 MCSE。
# Pairing: compare within the same MP and replicate (same random numbers); average within replicate,
# then take the mean and Monte Carlo standard error (MCSE) across replicates.

source("R/00_settings.R")
FACTORS <- c("S1", "S2", "S3", "S4")
dest <- file.path(OUT, "08_analysis"); dir.create(dest, showWarnings = FALSE)

# ---- 1. 读入 12 个配置 | read the 12 configurations ----
metrics_all <- list()
for (cid in CONFIGS$config_id) {
  f <- file.path(OUT, paste0("07_", cid), "replicate_metrics.csv.gz")
  if (!file.exists(f)) stop("先运行 07 | run step 07 first: CONFIG=", cid)
  x <- read.csv(gzfile(f)); x$config_id <- cid
  metrics_all[[cid]] <- x
}
all_metrics <- bind_rows(metrics_all)
cat("replicates per configuration:", max(all_metrics$rep), "| trajectories:", nrow(all_metrics), "\n")

# ---- 2. 绝对水平：每个配置 × 情景 × MP 的均值（补充表 7、9、10，图 5）----
# ---- 2. absolute levels per configuration × scenario × MP (Supplementary Tables 7, 9 and 10, Fig. 5) ----
num <- setdiff(names(all_metrics)[sapply(all_metrics, is.numeric)], c("rep", FACTORS))
absolute <- all_metrics |>
  pivot_longer(all_of(num), names_to = "metric", values_to = "value") |>
  group_by(config_id, case_id, mp, S1, S2, S3, S4, metric) |>
  # 未受配额约束的船队没有配额指标（NA），跳过 | fleets without a quota have no quota metrics (NA); skip them
  summarise(estimate = mean(value, na.rm = TRUE), n_rep = sum(!is.na(value)),
            mcse = sd(value, na.rm = TRUE) / sqrt(sum(!is.na(value))), .groups = "drop")

# ---- 3. 相对参考情景 C05（S1–S4 全关）的配对效应（补充图 2、3，补充表 15）----
# ---- 3. paired effects against C05 (Supplementary Figs 2 and 3, Supplementary Table 15) ----
outcomes <- c("mean_D", "D_overlap", "relative_catch", "mae_target", "outside_W_fraction")
paired <- list()
for (o in outcomes) {
  ref <- all_metrics[all_metrics$case_id == "C05", c("config_id", "mp", "rep", o)]
  names(ref)[4] <- "ref"
  z <- left_join(all_metrics[c("config_id", "case_id", "mp", "rep", FACTORS, o)], ref, by = c("config_id", "mp", "rep"))
  z$delta <- z[[o]] - z$ref
  paired[[o]] <- z |> group_by(config_id, case_id, mp, S1, S2, S3, S4) |>
    summarise(estimate = mean(delta), n_rep = n(), mcse = sd(delta) / sqrt(n()), .groups = "drop") |>
    mutate(metric = o)
}
paired <- bind_rows(paired)

# ---- 4. 边际效应：某因子开 vs 关，其他因子固定，在重复内对其他因子平均（图 4，补充表 8、14）----
# ---- 4. marginal effects: factor on vs off with the others fixed, averaged over the others within replicate (Fig. 4, Supplementary Tables 8 and 14) ----
marginal <- list()
for (f in FACTORS) for (o in c("mean_D", "D_overlap", "relative_catch", "mae_target")) {
  keys <- c("config_id", "mp", "rep", setdiff(FACTORS, f))
  base <- all_metrics[all_metrics[[f]] == 0, c(keys, o)]; names(base)[ncol(base)] <- "base"
  z <- inner_join(all_metrics[all_metrics[[f]] != 0, c(keys, f, o)], base, by = keys)
  z$delta <- z[[o]] - z$base
  marginal[[paste(f, o)]] <- z |>
    group_by(config_id, mp, rep, level = .data[[f]]) |> summarise(delta = mean(delta), .groups = "drop") |>   # 重复内 | within replicate
    group_by(config_id, mp, level) |>
    summarise(estimate = mean(delta), n_rep = n(), mcse = sd(delta) / sqrt(n()), .groups = "drop") |>
    mutate(contrast = paste0(f, ifelse(level > 0, "+", "-")), metric = o) |> select(-level)
}
marginal <- bind_rows(marginal)

# ---- 5. 估计方法得到真实移动信息时的变化（补充图 4）| change when the EM receives realized movement (Supplementary Fig. 4) ----
info_metrics <- c("mean_D", "relative_catch", "estimation_rmse_fixed", "estimation_bias_fixed")
info <- inner_join(metrics_all$movement_information[c("mp", "rep", FACTORS, info_metrics)],
                   metrics_all$reference[c("mp", "rep", FACTORS, info_metrics)],
                   by = c("mp", "rep", FACTORS), suffix = c("_observed", "_mean"))
information <- list()
for (o in info_metrics) {
  info$delta <- info[[paste0(o, "_observed")]] - info[[paste0(o, "_mean")]]
  information[[o]] <- info |> group_by(mp, rep, S3) |> summarise(delta = mean(delta), .groups = "drop") |>
    group_by(mp, S3) |> summarise(estimate = mean(delta), n_rep = n(), mcse = sd(delta) / sqrt(n()), .groups = "drop") |>
    mutate(metric = o)
}
information <- bind_rows(information)

# ---- 6. 保存（±1.96 MCSE）| save (±1.96 MCSE) ----
for (tab in c("absolute", "paired", "marginal", "information")) {
  x <- get(tab); x$lower <- x$estimate - 1.96 * x$mcse; x$upper <- x$estimate + 1.96 * x$mcse
  write.csv(x, file.path(dest, paste0(tab, ".csv")), row.names = FALSE)
}

# ---- 7. 正文中的数字 | numbers quoted in the Article ----
cat("\n== marginal effects, reference case: SB/SB_F=0 (percentage points) and catch (% of final historical catch) ==\n")
print(as.data.frame(marginal |> filter(config_id == "reference", metric %in% c("mean_D", "relative_catch")) |>
  mutate(value = round(100 * estimate, 2), mp = substr(mp, 1, 4)) |>
  select(metric, contrast, mp, value) |> pivot_wider(names_from = mp, values_from = value) |> arrange(metric, contrast)))

cat("\n== Fig. 5: MP-B with S1, S3 and S4 off (C01 under-reporting, C05 accurate, C09 over-reporting), percent ==\n")
print(as.data.frame(absolute |> filter(config_id == "reference", mp == "MP-B joint TAC", case_id %in% c("C01", "C05", "C09"),
                                       metric %in% c("mean_D", "relative_catch", "mae_target")) |>
  mutate(value = round(100 * estimate, 2)) |> select(case_id, metric, value) |> pivot_wider(names_from = metric, values_from = value)))

cat("\n== Supplementary Fig. 3, Table 15: combined scenario C23 vs C05, change in SB/SB_F=0 (percentage points) ==\n")
print(as.data.frame(paired |> filter(case_id == "C23", metric == "mean_D") |>
  mutate(value = round(100 * estimate, 2), mp = substr(mp, 1, 4)) |> select(config_id, mp, value) |>
  pivot_wider(names_from = mp, values_from = value)))

cat("\ntrajectories ever below the 20% limit reference point:", sum(all_metrics$P_ever_LRP),
    "of", nrow(all_metrics), "| lowest mean SB/SB_F=0:", round(min(all_metrics$mean_D), 4), "\n")

# ---- 8. 补充表 9：参考配置下时变移动（S3）关与开时的配额使用 | Supplementary Table 9: quota use with S3 off and on, reference case ----
# 第 2 节的均值（第 11–20 年与重复）再对 S1、S2、S4 的 12 种组合平均；未捕获的 TAC 是总配额中没有被捕获的百分比。
# The means of section 2 (years 11–20 and replicates) are averaged over the 12 combinations of S1, S2 and S4; the
# uncaught part of the TAC is the percentage of the total quota that is not caught.
t9 <- absolute |>
  filter(config_id == "reference", metric %in% c("quota_fully_caught_W", "quota_fully_caught_OW", "uncaught_quota_ratio", "relative_catch")) |>
  group_by(mp, S3, metric) |> summarise(value = 100 * mean(estimate), n_cells = n(), .groups = "drop")
stopifnot(all(t9$n_cells == 12))                       # 2 × 3 × 2 个组合 | 2 × 3 × 2 combinations
table9 <- t9 |> select(-n_cells) |> pivot_wider(names_from = metric, values_from = value) |> arrange(mp, S3) |>
  transmute(MP = substr(mp, 1, 4), `S3 Time-varying movement` = ifelse(S3 == 1, "On", "Off"),
            `Years in which fleet W takes its full quota (%)` = sprintf("%.1f", quota_fully_caught_W),
            `Years in which fleet OW takes its full quota (%)` = sprintf("%.1f", quota_fully_caught_OW),
            `Uncaught part of the TAC (%)` = sprintf("%.2f", uncaught_quota_ratio),
            `Catch (% of final historical catch)` = sprintf("%.1f", relative_catch))
write.csv(table9, file.path(dest, "supplementary_table_9.csv"), row.names = FALSE)
cat("\n== Supplementary Table 9: quota use with S3 off and on, reference case (percent) ==\n")
print(setNames(as.data.frame(table9), c("MP", "S3", "W full quota", "OW full quota", "uncaught TAC", "catch")), right = FALSE)

# ---- 9. 补充表 14：12 个配置中各因子对 SB/SB_F=0 的效应（百分点）----
# ---- 9. Supplementary Table 14: effect of each factor on SB/SB_F=0 in the 12 configurations (percentage points) ----
# 第 4 节的边际效应；四舍五入为 0 的值不带正负号 | marginal effects of section 4; values that round to zero are unsigned
config_labels <- c(reference = "Reference case", adult_M_0.24 = "Adult M = 0.24 yr⁻¹", adult_M_0.36 = "Adult M = 0.36 yr⁻¹",
                   adult_M_0.54 = "Adult M = 0.54 yr⁻¹", steepness_0.65 = "Steepness h = 0.65", steepness_0.95 = "Steepness h = 0.95",
                   connectivity_0.5 = "Connectivity × 0.5", connectivity_2 = "Connectivity × 2", selectivity_2.5 = "Selectivity a50 = 2.5",
                   selectivity_4.5 = "Selectivity a50 = 4.5", errors = "Errors in history and projection",
                   movement_information = "Movement information in estimation")
contrast_labels <- c("S1+" = "S1 Choice of applicable measures", "S2-" = "S2 Under-reporting", "S2+" = "S2 Over-reporting",
                     "S3+" = "S3 Time-varying movement", "S4+" = "S4 Catch outside the MP")
t14 <- marginal |> filter(metric == "mean_D") |>
  mutate(value = sprintf("%+.2f", 100 * estimate), value = ifelse(as.numeric(value) == 0, "0.00", value))
table14 <- t14 |>
  transmute(Configuration = factor(config_labels[config_id], config_labels), MP = substr(mp, 1, 4),
            contrast = factor(contrast_labels[contrast], contrast_labels), value) |>
  pivot_wider(names_from = contrast, values_from = value, names_sort = TRUE) |> arrange(Configuration, MP)
write.csv(table14, file.path(dest, "supplementary_table_14.csv"), row.names = FALSE)
cat("\n== Supplementary Table 14: effects on SB/SB_F=0 (percentage points) in the 12 configurations ==\n")
print(setNames(as.data.frame(table14), c("configuration", "MP", "S1", "S2 under", "S2 over", "S3", "S4")), right = FALSE)
# 表注中的两项检查 | the two checks stated in the table note
tested <- t14 |> filter((contrast %in% c("S2-", "S2+") & mp %in% c("MP-B joint TAC", "MP-C TAC with deduction")) | contrast %in% c("S3+", "S4+"))
cat("largest MCSE:", round(100 * max(t14$mcse), 3), "percentage points | reporting error under MP-B and MP-C, S3 and S4:",
    sum(abs(tested$estimate) > 1.96 * tested$mcse), "of", nrow(tested), "intervals of ±1.96 MCSE exclude zero\n")
cat("08 done:", dest, "\n")
