# 07_experiment.R — 一个模型配置的析因实验：24 个情景 × 4 个候选 MP × N_REP 个重复
# 07_experiment.R — factorial experiment for one model configuration: 24 scenarios × 4 candidate MPs × N_REP replicates
#
#   CONFIG=reference N_REP=20 Rscript R/07_experiment.R     # 快速查看 | quick look
#   CONFIG=reference N_REP=300 Rscript R/07_experiment.R    # 论文设置（约 6 分钟）| as in the Article (about 6 min)
# 12 个配置各运行一次（README 有并行命令）。重复 r 的随机数只取决于种子与 r。
# Run once per configuration (parallel command in the README). Replicate r depends only on the seed and r.

if (!exists("CONFIG")) CONFIG <- Sys.getenv("CONFIG", "reference")
if (!exists("N_REP"))  N_REP  <- as.integer(Sys.getenv("N_REP", "20"))
RUN_CELLS <- NULL
source("R/00_settings.R")
RUN_CELLS <- design$case_id                           # 全部 24 个情景 | all 24 scenarios
RUN_MPS <- c("MP-A WCPFC TAC", "MP-B joint TAC", "MP-C TAC with deduction", "MP-D overlap effort cap")
FIXED_SHARES <- FALSE

# ---- 1. 依次运行 03–06 | run steps 03–06 in order ----
t0 <- Sys.time()
source("R/03_parameters.R")                           # 参数 | parameters
source("R/04_initialise.R")                           # 平衡与校准 | equilibrium and calibration
source("R/05_history.R")                              # 随机数与历史 | random numbers and history
source("R/06_closed_loop.R")                          # 闭环预测 | closed-loop projection
cat("elapsed:", format(Sys.time() - t0), "\n")

# ---- 2. 保存 | save ----
dest <- file.path(OUT, paste0("07_", CONFIG))
dir.create(dest, showWarnings = FALSE)
con <- gzfile(file.path(dest, "replicate_metrics.csv.gz"), "wt")
write.csv(metrics, con, row.names = FALSE, na = ""); close(con)
write.csv(trajectories, file.path(dest, "first_replicate_trajectories.csv"), row.names = FALSE, na = "")
write.csv(diagnostics, file.path(dest, "diagnostics.csv"), row.names = FALSE, na = "")
cat("07 done:", dest, "\n")
