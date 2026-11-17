# 03_parameters.R — 操作模型、估计方法与 HCR 的全部参数（表 1、补充表 9–10）
# 03_parameters.R — all parameters of the operating model, estimation method and HCR (Table 1, Supplementary Tables 9–10)

source("R/00_settings.R")

# ---- 1. 生活史：12 个年龄组，末组 12+（SPC 2024 诊断模型）| life history: 12 age classes, last is 12+ (SPC 2024 diagnostic model) ----
bio <- fromJSON(file.path(DATA, "biology_at_age.json"))
mort <- read.csv(file.path(DATA, "natural_mortality_steepness.csv"))
mort <- mort[mort$biology == cfg$biology, ]             # 本配置的 M 与 h | M and h of this configuration
age <- as.integer(bio$age)
A   <- length(age)
w   <- bio$weight_t                                     # 年龄别平均体重（吨）| mean weight at age (t)
mat <- bio$reproductive_potential                       # 性比 × 成熟度 | sex ratio × maturity
M   <- mort$M                                           # 年龄别自然死亡（年⁻¹）| natural mortality at age (yr⁻¹)
h   <- unique(mort$h)                                   # Beverton–Holt 陡度 | steepness
R0  <- 1e6                                              # 未捕捞补充量（数值尺度）| unfished recruitment (numerical scale)
stopifnot(A == 12, length(M) == A, length(h) == 1)

# ---- 2. 选择性：所有船队共用 logistic，中点 3.5 岁（敏感性 2.5、4.5）| selectivity: one logistic for all fleets, midpoint 3.5 (2.5, 4.5) ----
midpoint <- c(common = 3.5, younger = 2.5, older = 4.5)[[cfg$selectivity]]
sel <- plogis((age - midpoint) / .6)

# ---- 3. 三个区域与移动 | three regions and movement ----
# R1 仅 WCPFC 公约区，R2 重疊区，R3 仅 IATTC 公约区（东太）| R1 WCPFC only, R2 overlap area, R3 EPO only
recruit_share <- c(.8, .1, .1)                          # 补充量的区域分配 | regional split of recruitment
# 年移动概率：行 = 出发区，列 = 到达区 | annual movement probabilities: rows = origin, columns = destination
movement_base <- matrix(c(.90, .10, 0,
                          .15, .75, .10,
                          0,   .15, .85), 3, byrow = TRUE)
k_conn <- cfg$connectivity                              # 连通度倍数（0.5、1、2）| connectivity scale (0.5, 1, 2)
movement <- diag(3) + k_conn * (movement_base - diag(3))  # P(k) = I + k (P − I)，方程 4 | equation 4
# S3：环境状态 −1/0/+1 使 R1→R2、R2→R3 的概率各变 ±0.04、±0.03
# S3: state −1/0/+1 shifts the R1→R2 and R2→R3 probabilities by ±0.04 and ±0.03
movement_shift <- k_conn * c(.04, .03)
env_transition <- rbind(c(.4, .5, .1), c(.25, .5, .25), c(.1, .5, .4))   # 行 = 上一年状态 | rows = previous state
estimation_movement <- cfg$estimation_movement          # mean：估计方法用平均移动 | mean: the EM uses mean movement

# ---- 4. 四支船队 | four fleets ----
# W：R1，WCPFC 措施；OW、OI：R2，分别适用 WCPFC、IATTC 措施；E：R3
# W: R1 under WCPFC measures; OW and OI: R2 under WCPFC and IATTC measures; E: R3
fleets <- c("W", "OW", "OI", "E")
fleet_area <- c(1L, 2L, 2L, 3L)
effort_end <- c(W = 1, OW = .2, OI = .1, E = .3)        # 历史末年的相对容量 | relative capacity at the end of the history
effort_creep <- 1.01                                    # 可捕系数每年增长 1% | catchability growth 1% per year
epo_growth <- .03                                       # S4：E 容量每年再增 3% | S4: E capacity grows a further 3% per year

# ---- 5. 时间 | time ----
n_hist <- 30L; n_proj <- 20L                            # 合成历史 30 年，预测 20 年 | 30-year history, 20-year projection
history_years <- 1995:2024                              # 历史年的日历标签 | calendar labels of the history
evaluation_years <- 11:20                               # 评估期 = 预测第 11–20 年 | evaluation period
cycle <- 3L                                             # 每 3 年一次决策（CMM 2025-01）| decision every 3 years

# ---- 6. HCR（CMM 2025-01 附件 I）| harvest control rule (CMM 2025-01 Annex I) ----
# 控制点（相对状态，渔获乘数）| control points (relative status, catch scalar)
hcr_points <- rbind(c(.37, .20), c(.94, 1.09), c(1.29, 1.09), c(1.59, 1.29))
hcr_curvature <- 10                                     # 0.37–0.94 之间的指数段 | exponential section between 0.37 and 0.94
change_bounds <- c(.95, 1.10)                           # 每次决策 −5%/+10% | −5%/+10% per decision
state_reference_years <- 2017:2019                      # 相对状态的分母 | denominator of relative status
catch_reference_years <- 2020:2022                      # 渔获乘数的基数 | base of the catch scalar
target_relative <- .96                                  # 目标 = 0.96 × 2017–2019 耗竭 | target = 0.96 × 2017–2019 depletion
unfished_window <- 10L                                  # 未捕捞参照：前 10 年平均 | unfished reference: preceding 10-year mean
limit_reference <- .20                                  # 限制参考点 | limit reference point
diagnostic_threshold <- .42                             # 诊断阈值 | diagnostic threshold

# ---- 7. 初始化与估计方法 | initialisation and estimation method ----
start_depletion <- .45                                  # 历史末年的耗竭 | depletion at the end of the history
illustrative_target <- .56                              # 定参考船队份额的平衡耗竭 | equilibrium depletion for reference fleet shares
index_gain <- .5                                        # 估计方法的指数修正权重 | index weight in the estimation method
em_Fmax <- 6                                            # 估计方法的 F 上限 | F bound in the estimation method

# ---- 8. 实施因子 | implementation factors ----
misreporting <- .25                                     # S2：中心 ±25%，每年乘 U(0.5, 1.5) | centre ±25%, times U(0.5, 1.5) each year
n_segments <- 60L                                       # 重疊区 60 个决策单元 | 60 decision units in the overlap area
min_measure_years <- 3L                                 # 至少 3 年才能更换（IATTC C-12-11）| at least 3 years between changes
switch_cost <- .05; switch_intercept <- -2; switch_slope <- 6   # 更换措施的 logit | logit of changing measures
switch_strength <- .3                                   # 概率上限 | probability ceiling
overlap_cap_fraction <- .9448275862068966               # MP-D 重疊区努力量上限 / 容量 | MP-D overlap effort cap / capacity

# ---- 9. 误差（参考配置为零）| errors (zero in the reference case) ----
# 对数正态乘子 exp(σz − σ²/2)，均值为 1 | lognormal multipliers with mean one
if (cfg$errors == "illustrative") {
  sigma_R <- .4; sigma_catch <- .05; sigma_index <- .15; sigma_impl <- .10
} else {
  sigma_R <- 0;  sigma_catch <- 0;   sigma_index <- 0;   sigma_impl <- 0
}

# ---- 10. 四个候选 MP | four candidate MPs ----
mps <- c("MP-A WCPFC TAC", "MP-B joint TAC", "MP-C TAC with deduction", "MP-D overlap effort cap")

# ---- 11. 查看 | look ----
print(data.frame(age, weight_t = w, maturity = round(mat, 4), M, selectivity = round(sel, 4)))
cat("steepness h =", h, "| adult mean M (ages 4–12) =", round(mean(M[4:12]), 4), "\n")
print(movement)
cat("03 done: CONFIG =", CONFIG, "\n")
