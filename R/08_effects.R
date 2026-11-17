# 08_effects.R — 配对效应、边际效应与正文中的数字（图 4–6、补充表 4–8）
# 08_effects.R — paired effects, marginal effects and the numbers quoted in the Article (Figs 4–6,
#                Supplementary Tables 4–8)
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

# ---- 2. 绝对水平：每个配置 × 情景 × MP 的均值（补充表 4、图 5）| absolute levels per configuration × scenario × MP (Supplementary Table 4, Fig. 5) ----
num <- setdiff(names(all_metrics)[sapply(all_metrics, is.numeric)], c("rep", FACTORS))
absolute <- all_metrics |>
  pivot_longer(all_of(num), names_to = "metric", values_to = "value") |>
  group_by(config_id, case_id, mp, S1, S2, S3, S4, metric) |>
  # 未受配额约束的船队没有配额指标（NA），跳过 | fleets without a quota have no quota metrics (NA); skip them
  summarise(estimate = mean(value, na.rm = TRUE), n_rep = sum(!is.na(value)),
            mcse = sd(value, na.rm = TRUE) / sqrt(sum(!is.na(value))), .groups = "drop")

# ---- 3. 相对参考情景 C05（S1–S4 全关）的配对效应（补充图 2、图 6）| paired effects against C05 (Supplementary Fig. 2, Fig. 6) ----
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

# ---- 4. 边际效应：某因子开 vs 关，其他因子固定，在重复内对其他因子平均（图 4、补充表 5）----
# ---- 4. marginal effects: factor on vs off with the others fixed, averaged over the others within replicate (Fig. 4, Supplementary Table 5) ----
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

# ---- 5. 估计方法得到真实移动信息时的变化（补充图 3）| change when the EM receives realized movement (Supplementary Fig. 3) ----
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
cat("\n== marginal effects, reference case: depletion (percentage points) and catch (% of final historical catch) ==\n")
print(as.data.frame(marginal |> filter(config_id == "reference", metric %in% c("mean_D", "relative_catch")) |>
  mutate(value = round(100 * estimate, 2), mp = substr(mp, 1, 4)) |>
  select(metric, contrast, mp, value) |> pivot_wider(names_from = mp, values_from = value) |> arrange(metric, contrast)))

cat("\n== Fig. 5: MP-B with misreporting only (C01 under, C05 unbiased, C09 over), percent ==\n")
print(as.data.frame(absolute |> filter(config_id == "reference", mp == "MP-B joint TAC", case_id %in% c("C01", "C05", "C09"),
                                       metric %in% c("mean_D", "relative_catch", "mae_target")) |>
  mutate(value = round(100 * estimate, 2)) |> select(case_id, metric, value) |> pivot_wider(names_from = metric, values_from = value)))

cat("\n== Fig. 6: combined scenario C23 vs C05, change in depletion (percentage points) ==\n")
print(as.data.frame(paired |> filter(case_id == "C23", metric == "mean_D") |>
  mutate(value = round(100 * estimate, 2), mp = substr(mp, 1, 4)) |> select(config_id, mp, value) |>
  pivot_wider(names_from = mp, values_from = value)))

cat("\ntrajectories ever below the 20% limit reference point:", sum(all_metrics$P_ever_LRP),
    "of", nrow(all_metrics), "| lowest mean depletion:", round(min(all_metrics$mean_D), 4), "\n")
cat("08 done:", dest, "\n")
