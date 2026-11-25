# 00_settings.R — 路径、程序包、随机种子、模型配置与析因设计
# 00_settings.R — paths, packages, random seed, model configurations and factorial design
# 在仓库根目录运行；其余脚本都会先 source 本文件。
# Run from the repository root; every other script sources this file first.

# ---- 1. 程序包 | packages ----
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(jsonlite)
})

# ---- 2. 路径 | paths ----
DATA <- "data"                                       # 输入数据 | input data
RAW  <- file.path("data", "raw")                     # 公开延绳钓原始数据 | raw public longline data
OUT  <- "output"                                     # 结果 | results
FIG  <- "figures"                                    # 图 | figures
dir.create(OUT, showWarnings = FALSE); dir.create(FIG, showWarnings = FALSE)

# ---- 3. 随机种子 | random seed ----
# 同一重复在各情景、各候选 MP 之间共用随机数，因此可以配对比较。
# Replicate r shares its random numbers across scenarios and candidate MPs, so outcomes are paired.
SEED <- 201306L

# ---- 4. 12 个模型配置 | the 12 model configurations ----
# reference = 参考配置；其余各改一项（表 2、补充图 3）| reference case; each other configuration changes one item (Table 2, Supplementary Fig. 3)
CONFIGS <- tribble(
  ~config_id,             ~biology,          ~errors,        ~connectivity, ~selectivity, ~estimation_movement,
  "reference",            "reference",       "none",         1,             "common",     "mean",
  "adult_M_0.24",         "adult_M_0.24",    "none",         1,             "common",     "mean",
  "adult_M_0.36",         "adult_M_0.36",    "none",         1,             "common",     "mean",
  "adult_M_0.54",         "adult_M_0.54",    "none",         1,             "common",     "mean",
  "steepness_0.65",       "steepness_0.65",  "none",         1,             "common",     "mean",
  "steepness_0.95",       "steepness_0.95",  "none",         1,             "common",     "mean",
  "connectivity_0.5",     "reference",       "none",         0.5,           "common",     "mean",
  "connectivity_2",       "reference",       "none",         2,             "common",     "mean",
  "selectivity_2.5",      "reference",       "none",         1,             "younger",    "mean",
  "selectivity_4.5",      "reference",       "none",         1,             "older",      "mean",
  "errors",               "reference",       "illustrative", 1,             "common",     "mean",
  "movement_information", "reference",       "none",         1,             "common",     "observed")

# ---- 5. 本次运行的选择（可在 source 之前赋值）| choices for this run (assign before sourcing to override) ----
if (!exists("CONFIG")) CONFIG <- Sys.getenv("CONFIG", "reference")
if (!exists("N_REP"))  N_REP  <- as.integer(Sys.getenv("N_REP", "2"))   # 论文用 300 | 300 in the Article
stopifnot(CONFIG %in% CONFIGS$config_id, N_REP >= 1)
cfg <- CONFIGS[CONFIGS$config_id == CONFIG, ]

# ---- 6. 析因设计：2 × 3 × 2 × 2 = 24 个情景 | factorial design: 2 × 3 × 2 × 2 = 24 scenarios ----
# S1 选择适用措施（0/1）；S2 重疊区渔获的上报（−1 少报、0 准确、+1 多报）；S3 时变移动（0/1）；
# S4 MP 范围外渔获，即 E 船队容量增长（0/1）
# S1 choice of applicable measures (0/1); S2 reporting of overlap catch (−1 under-reporting, 0 accurate, +1 over-reporting);
# S3 time-varying movement (0/1); S4 catch outside the MP, as growth of fleet-E capacity (0/1)
design <- expand.grid(S4 = 0:1, S3 = 0:1, S2 = -1:1, S1 = 0:1)
design <- data.frame(case_id = sprintf("C%02d", 1:24), S1 = design$S1, S2 = design$S2, S3 = design$S3, S4 = design$S4)
# C05 = 参考情景（S1–S4 全关）；C23 = 合并情景（S1 开、多报、S3 开）
# C05 = reference scenario (S1–S4 off); C23 = combined scenario (S1 on, over-reporting, S3 on)

cat("00 settings: CONFIG =", CONFIG, "| N_REP =", N_REP, "| seed =", SEED, "\n")
