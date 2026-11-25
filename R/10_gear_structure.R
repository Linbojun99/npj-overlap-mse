# 10_gear_structure.R — 附加分析：W 船队拆成延绳钓（LL）与曳绳钓（TR）；论文的图表不用这一步的结果
# 10_gear_structure.R — additional analysis: the W fleet split into longline (LL) and troll (TR); no display of the
#                       Article uses the results of this step
#
# 六种结构 | six structures:
#   shared_2pct / shared_5pct     TR 占历史 W+OW 渔获 2% / 5%，与 LL 选择性相同 | TR 2% / 5% of historical W+OW catch, LL selectivity
#   troll11_2pct / troll11_5pct   同上，TR 用 SPC 诊断模型渔业 11 的选择性 | as above, selectivity of fishery 11 (SPC diagnostic model)
#   troll12_2pct / troll12_5pct   同上，渔业 12 | as above, fishery 12
# 没有 TR 的四船队结构就是 07 的参考配置 | the four-fleet structure without troll is step 07's reference case
#   GEAR_CONFIG=troll11_5pct N_REP=20 Rscript R/10_gear_structure.R
#   GEAR_CONFIG=troll11_5pct N_REP=300 Rscript R/10_gear_structure.R    # 300 个重复（约 30 分钟）| 300 replicates (about 30 min)

# ---- 0. 本次运行的选择 | choices for this run ----
if (!exists("GEAR_CONFIG")) GEAR_CONFIG <- Sys.getenv("GEAR_CONFIG", "troll11_5pct")
if (!exists("N_REP"))       N_REP       <- as.integer(Sys.getenv("N_REP", "20"))
CONFIG <- "reference"                                  # 参考配置的生物学 | reference-case biology
t0 <- Sys.time()

# ---- 1. 参数（03）、平衡（04）、随机数（05）| parameters (03), equilibrium (04), random numbers (05) ----
source("R/03_parameters.R")
source("R/04_initialise.R")
source("R/05_history.R")                               # draws 用于下面的历史与预测 | draws are reused below

# ---- 2. 渔具结构 | gear structure ----
# TR 份额 = TR /（W + OW）历史渔获重量；2%、5% 取自 SC22-SA-IP15 表 1–2 南太平洋 LL+TR 比例（取整）
# TR share = TR / (W + OW) historical catch weight; 2% and 5% round the southern LL+TR ratios of SC22-SA-IP15 Tables 1–2
GEAR_CONFIGS <- tribble(
  ~config_id,     ~share, ~type,
  "shared_2pct",  .02,    "shared",
  "shared_5pct",  .05,    "shared",
  "troll11_2pct", .02,    "fishery11",
  "troll11_5pct", .05,    "fishery11",
  "troll12_2pct", .02,    "fishery12",
  "troll12_5pct", .05,    "fishery12")
stopifnot(GEAR_CONFIG %in% GEAR_CONFIGS$config_id)
gear_share <- GEAR_CONFIGS$share[GEAR_CONFIGS$config_id == GEAR_CONFIG]
gear_type  <- GEAR_CONFIGS$type[GEAR_CONFIGS$config_id == GEAR_CONFIG]
dest <- file.path(OUT, paste0("10_", GEAR_CONFIG)); dir.create(dest, showWarnings = FALSE)

# 五种渔具：4 个延绳钓段 + R1 的曳绳钓 | five gears: the 4 longline segments + troll in R1
gear_names <- c("W_LL", "OW", "OI", "E", "TR")
gear_area  <- c(1L, 2L, 2L, 3L, 1L)
# TR 选择性：shared = 与 LL 相同；fishery11/12 = SPC 诊断模型的拟合选择性 | TR selectivity: LL curve, or fitted fisheries 11/12
tr_sel <- sel
if (gear_type != "shared") tr_sel <- read.csv(file.path(DATA, "troll_selectivity.csv"))[[gear_type]]
gear_sel <- rbind(W_LL = sel, OW = sel, OI = sel, E = sel, TR = tr_sel)   # 5 × 12
S_R1 <- gear_sel[c("W_LL", "TR"), ]                    # R1 的两种渔具 | the two gears of R1
same_sel <- all(S_R1 == matrix(S_R1[1, ], 2, A, byrow = TRUE))   # 选择性相同时联合求解是一维的 | identical rows: one-dimensional solve
TR_Fmax <- em_Fmax                                     # TR 的 F 上限（数值界，6 年⁻¹）| TR F ceiling (numerical bound, 6 yr⁻¹)
TR_Fcap <- TR_Fmax
cat("\ngear structure:", GEAR_CONFIG, "| TR share =", gear_share, "| TR selectivity =", gear_type, "\n")
print(data.frame(age, LL = round(sel, 4), TR = round(tr_sel, 4)), row.names = FALSE)

# ---- 3. 校准：LL 可捕系数 q 与 TR/LL 的 F 比 | calibration: LL catchability q and the TR/LL F ratio ----
# 两个条件同时满足 | two conditions hold together:
#   (a) 固定努力量下的平衡耗竭 = 0.45（同 04）| equilibrium depletion at fixed effort = 0.45 (as in 04)
#   (b) 平衡时 TR 渔获 /（W + OW）渔获 = TR 份额 | at equilibrium TR catch / (W + OW) catch = TR share
# 平衡状态用每补充量的年龄–区域数量与 Beverton–Holt 解析求出（不再跑 400 年）。
# The equilibrium is solved analytically from age–region numbers per recruit and Beverton–Holt (no 400-year loop).
# 平衡函数：给定五种渔具的 F，求平衡种群并走一年，返回产卵潜力与渔获（uniroot 需要函数）
# equilibrium function: for the F of five gears, solve the equilibrium and step one year (uniroot needs a function)
gear_equilibrium <- function(Fg) {
  Z <- rbind(Fg[1] * sel + Fg[5] * tr_sel, (Fg[2] + Fg[3]) * sel, Fg[4] * sel) + Mmat   # R1 两种渔具的死亡率相加 | R1 adds both gears
  P <- t(movement)
  nr <- matrix(0, 3, A); nr[, 1] <- recruit_share   # 每个补充量的年初数量 | start-of-year numbers per recruit
  for (a in 2:(A - 1)) nr[, a] <- exp(-Z[, a - 1]) * as.numeric(P %*% nr[, a - 1])   # 移动后存活到下一龄 | move, then survive
  # 12+ 组：n = s11 P n11 + s12 P n，解线性方程 | plus group: solve n = s11 P n11 + s12 P n
  nr[, A] <- solve(diag(3) - diag(exp(-Z[, A])) %*% P, exp(-Z[, A - 1]) * as.numeric(P %*% nr[, A - 1]))
  phi <- sum(as.numeric(nr %*% (mat * w)))           # 每补充量产卵潜力 | spawning potential per recruit
  Rq <- max(0, (4 * h * R0 * phi - SB0 * (1 - h)) / (phi * (5 * h - 1)))   # Beverton–Holt 平衡补充 | equilibrium recruitment
  N <- nr * Rq
  Nm <- t(movement) %*% N                            # 再走一年：移动 → 捕捞 | one more year: move → fish
  numbers <- gear_sel * (Nm * (-expm1(-Z)) / Z)[gear_area, ] * Fg   # 各渔具渔获尾数 5 × 12 | catch numbers by gear
  cg <- as.numeric(numbers %*% w)                    # 渔获重量 | catch weight
  list(state = N, SB = as.numeric(Nm %*% (mat * w)), catch_gear = cg,
       catch = c(W = cg[1] + cg[5], OW = cg[2], OI = cg[3], E = cg[4]))   # 母船队 W = LL + TR | parent W = LL + TR
}
# 给定 TR/LL 的 F 比：解 q 使耗竭 = 0.45，再算 TR 份额（外层 uniroot 需要函数）
# for a given TR/LL F ratio: solve q for depletion 0.45, then the TR share (the outer uniroot needs a function)
calibrate_ratio <- function(ratio) {
  qq <- uniroot(function(qx) { F4 <- qx * effort_end; sum(gear_equilibrium(c(F4, F4[1] * ratio))$SB) / SB0 - start_depletion },
                c(.00001, 3), tol = 1e-12)$root
  F4 <- qq * effort_end; o <- gear_equilibrium(c(F4, F4[1] * ratio))
  list(q = qq, o = o, share = o$catch_gear[5] / sum(o$catch[c("W", "OW")]))
}
# 根区间是 F 比；上界从 1 起翻倍到包住目标 | the bracket is an F ratio; double the upper end until it brackets the target
upper <- 1
while (calibrate_ratio(upper)$share < gear_share && upper < 128) upper <- upper * 2
tr_ratio <- uniroot(function(rr) calibrate_ratio(rr)$share - gear_share, c(0, upper), tol = 1e-12)$root
fit <- calibrate_ratio(tr_ratio)
q_end <- fit$q                                       # 历史末年 LL 可捕系数 | LL catchability at the end of the history
q_initial <- q_end / effort_creep^(n_hist - 1)       # 历史第 1 年 | first history year
N_start <- fit$o$state                               # 历史第 1 年的种群 | population at the start of the history
tr_alloc <- unname(fit$o$catch_gear[5] / fit$o$catch["W"])   # W 母配额中 TR 的固定比例 | fixed TR fraction of the parent W quota
realised_share <- unname(fit$share)
# 参考船队份额：平衡耗竭 0.56（同 04）| reference shares: equilibrium depletion 0.56 (as in 04)
scale_gear <- uniroot(function(sx) { F4 <- q_end * effort_end * sx; sum(gear_equilibrium(c(F4, F4[1] * tr_ratio))$SB) / SB0 - illustrative_target },
                      c(.001, 3), tol = 1e-10)$root
F4 <- q_end * effort_end * scale_gear; eq_gear <- gear_equilibrium(c(F4, F4[1] * tr_ratio))
C_target <- sum(eq_gear$catch)
share_ref <- eq_gear$catch / C_target                # W（含 TR）、OW、OI、E | W (incl. TR), OW, OI, E
managed_share <- sum(share_ref[c("W", "OW")])
stopifnot(abs(sum(fit$o$SB) / SB0 - start_depletion) < 1e-8, abs(realised_share - gear_share) < 1e-9)
calibration <- data.frame(
  quantity = c("q_end", "q_initial", "TR_LL_F_ratio", "TR_fraction_of_W_quota", "TR_share", "C_target", paste0("share_ref_", fleets)),
  value = c(q_end, q_initial, tr_ratio, tr_alloc, realised_share, C_target, share_ref))
print(calibration, row.names = FALSE)
write.csv(calibration, file.path(dest, "calibration.csv"), row.names = FALSE)

# ---- 4. 30 年历史：同样的随机数，种群、观测与估计按五种渔具计算 | 30-year history: same draws, five gears ----
for (r in seq_len(N_REP)) {
  d <- draws[[r]]
  N <- N_start; est <- N_start; U <- N0              # 真实、估计、未捕捞配对种群 | true, estimated, paired unfished
  U_history <- rep(SB0, unfished_window)
  catch_hist <- matrix(0, n_hist, 4, dimnames = list(NULL, fleets))
  est_D <- true_D <- numeric(n_hist)
  for (yr in seq_len(n_hist)) {
    rmult <- exp(sigma_R * d$recruitment_z[yr] - .5 * sigma_R^2)    # 补充乘子（= 1）| recruitment multiplier (= 1)
    # -- 真实种群：TR 的 F = W 延绳钓 F × 校准的 F 比 | true population: TR F = W longline F × calibrated ratio --
    q <- q_initial * effort_creep^(yr - 1)
    Fy <- q * history_effort[yr, ]
    Fg <- c(Fy, TR = Fy[[1]] * tr_ratio)             # W_LL, OW, OI, E, TR
    Nm <- t(movement) %*% N
    SB <- as.numeric(Nm %*% (mat * w)); VB <- as.numeric(Nm %*% (sel * w))   # 指数仍按 LL 选择性 | the index keeps LL selectivity
    Z <- rbind(Fg[1] * sel + Fg[5] * tr_sel, (Fg[2] + Fg[3]) * sel, Fg[4] * sel) + Mmat
    numbers <- gear_sel * (Nm * (-expm1(-Z)) / Z)[gear_area, ] * Fg   # 各渔具渔获尾数 | catch numbers by gear
    catch_gear <- setNames(as.numeric(numbers %*% w), gear_names)
    catch <- c(W = catch_gear[["W_LL"]] + catch_gear[["TR"]], catch_gear[c("OW", "OI", "E")])   # 母船队 | parent fleets
    rec <- 4 * h * R0 * sum(SB) / (SB0 * (1 - h) + sum(SB) * (5 * h - 1)) * rmult
    surv <- Nm * exp(-Z)
    N <- matrix(0, 3, A); N[, 2:A] <- surv[, 1:(A - 1)]; N[, A] <- N[, A] + surv[, A]; N[, 1] <- rec * recruit_share
    # -- 未捕捞配对种群（同 05）| paired unfished population (as in 05) --
    Um <- t(movement) %*% U
    SBu <- as.numeric(Um %*% (mat * w))
    Zu <- outer(c(0, 0, 0), sel) + Mmat
    recu <- 4 * h * R0 * sum(SBu) / (SB0 * (1 - h) + sum(SBu) * (5 * h - 1)) * rmult
    survu <- Um * exp(-Zu)
    U <- matrix(0, 3, A); U[, 2:A] <- survu[, 1:(A - 1)]; U[, A] <- U[, A] + survu[, A]; U[, 1] <- recu * recruit_share
    # -- 观测：TR 与母船队 W 共用观测误差乘子（不新增随机数）| observations: TR shares the W multiplier (no new random numbers) --
    obs_catch <- catch * exp(sigma_catch * d$catch_z[yr, ] - .5 * sigma_catch^2)
    obs_TR <- catch_gear[["TR"]] * exp(sigma_catch * d$catch_z[yr, 1] - .5 * sigma_catch^2)
    gear_obs <- c(W_LL = obs_catch[["W"]] - obs_TR, obs_catch[c("OW", "OI", "E")], TR = obs_TR)   # W 拆回 LL 与 TR | W split back into LL and TR
    index <- VB * exp(sigma_index * d$index_z[yr, ] - .5 * sigma_index^2)
    area_catch <- c(sum(obs_catch[fleet_area == 1]), sum(obs_catch[fleet_area == 2]), sum(obs_catch[fleet_area == 3]))
    # -- 估计方法：知道两种渔具的选择性；R1 用分渔具渔获反推两个 F | EM: knows both selectivities; R1 back-calculates two F from gear catches --
    Nm_e <- t(movement) %*% est
    VB_pred <- pmax(as.numeric(Nm_e %*% (sel * w)), 1e-12)
    filtered <- Nm_e * exp(index_gain * log(pmax(index, 1e-12) / VB_pred))
    SB_e <- as.numeric(filtered %*% (mat * w))
    F_e <- setNames(numeric(5), gear_names)
    # R1 联合 Baranov：每种渔具要么达到渔获，要么到 F 上限
    # R1 joint Baranov: each gear either reaches its catch or its F ceiling
    Na <- filtered[1, ]; Tj <- gear_obs[c("W_LL", "TR")]; Fcap <- c(em_Fmax, TR_Fcap)
    if (sum(Na * w) == 0 || all(Tj == 0) || all(Fcap == 0)) {
      F_e[c(1, 5)] <- Fcap * 0
    } else if (same_sel) {
      # 选择性相同：两种渔具的 F 之和是一维问题，与 05 的求解相同 | same selectivity: total F is one-dimensional, as in 05
      req <- function(x) { Zx <- M + x * sel; pmin(Fcap, Tj / sum(sel / Zx * (-expm1(-Zx)) * Na * w)) }
      tot <- if (sum(req(sum(Fcap))) - sum(Fcap) >= -1e-12) sum(Fcap) else
        uniroot(function(x) sum(req(x)) - x, c(0, sum(Fcap)), tol = 1e-11)$root
      F_e[c(1, 5)] <- req(tot)
    } else {
      # 选择性不同：阻尼牛顿法同时解两个 F（最多 60 步，每步最多减半 25 次）
      # different selectivities: damped Newton for both F together (at most 60 steps, at most 25 halvings each)
      Nw <- Na * w
      res <- function(x) {                           # 牛顿法与线搜索需要残差函数 | Newton and its line search need the residual
        Zx <- M + as.numeric(x %*% S_R1)
        gx <- as.numeric(S_R1 %*% (Nw * (-expm1(-Zx)) / Zx))       # 每单位 F 的渔获 | catch per unit F
        x - pmin(Fcap, ifelse(gx > 0, Tj / gx, ifelse(Tj > 0, Inf, 0)))
      }
      g0 <- as.numeric(S_R1 %*% (Nw * (-expm1(-M)) / M))           # F = 0 处的起点 | starting point at F = 0
      F2 <- pmin(Fcap, ifelse(g0 > 0, Tj / g0, Inf)); F2[Tj == 0] <- 0
      for (it in 1:60) {
        r0 <- res(F2)
        if (max(abs(r0)) < 1e-12) break                             # 收敛 | converged
        Zx <- M + as.numeric(F2 %*% S_R1); gx <- as.numeric(S_R1 %*% (Nw * (-expm1(-Zx)) / Zx))
        free <- ifelse(gx > 0, Tj / gx, ifelse(Tj > 0, Inf, 0)) < Fcap & Tj > 0 & is.finite(Tj)   # 未到上限的渔具 | gears below their ceiling
        dg <- (S_R1 * rep(Nw * ((Zx + 1) * exp(-Zx) - 1) / Zx^2, each = 2)) %*% t(S_R1)       # dg/dF
        J <- diag(2)                                                # 雅可比矩阵 | Jacobian
        if (any(free)) J[free, ] <- J[free, , drop = FALSE] + dg[free, , drop = FALSE] * (Tj[free] / gx[free]^2)
        step <- as.numeric(solve(J, r0))
        hv <- 0                                                     # 步长减半直到残差下降 | halve the step until the residual falls
        while (max(abs(res(pmax(0, pmin(Fcap, F2 - step * 2^(-hv)))))) >= max(abs(r0))) {
          hv <- hv + 1; if (hv > 25) stop("Joint gear solver line search failed")
        }
        F2 <- pmax(0, pmin(Fcap, F2 - step * 2^(-hv)))
        if (it == 60) stop("Joint gear solver did not converge")
      }
      F_e[c(1, 5)] <- F2
    }
    for (a in 2:3) {                                 # R2、R3 与 05 相同（R2 放 OW 位置，R3 放 E 位置）| R2, R3 as in 05 (R2 in the OW slot, R3 in the E slot)
      jj <- c(2, 4)[a - 1]
      Na <- filtered[a, ]; Ta <- area_catch[a]
      if (sum(Na * w) == 0 || Ta == 0) next
      req <- function(x) { Zx <- M + x * sel; pmin(em_Fmax, Ta / sum(sel / Zx * (-expm1(-Zx)) * Na * w)) }
      tot <- if (sum(req(em_Fmax)) - em_Fmax >= -1e-12) em_Fmax else
        uniroot(function(x) sum(req(x)) - x, c(0, em_Fmax), tol = 1e-11)$root
      F_e[jj] <- req(tot)
    }
    Ze <- rbind(F_e[1] * sel + F_e[5] * tr_sel, (F_e[2] + F_e[3]) * sel, F_e[4] * sel) + Mmat
    rec_e <- 4 * h * R0 * sum(SB_e) / (SB0 * (1 - h) + sum(SB_e) * (5 * h - 1))
    surve <- filtered * exp(-Ze)
    est <- matrix(0, 3, A); est[, 2:A] <- surve[, 1:(A - 1)]; est[, A] <- est[, A] + surve[, A]; est[, 1] <- rec_e * recruit_share
    # -- 记录 | record --
    true_D[yr] <- sum(SB) / mean(tail(U_history, unfished_window))
    U_history <- c(U_history, sum(SBu))
    est_D[yr] <- sum(SB_e) / SB0
    catch_hist[yr, ] <- obs_catch                    # 母船队上报渔获（W 含 TR）| parent reported catch (W includes TR)
  }
  si <- match(state_reference_years, history_years); ci <- match(catch_reference_years, history_years)
  init[[r]] <- list(
    N = N, U = U, est = est, U_history = tail(U_history, unfished_window),
    state_reference = mean(est_D[si]), true_state_reference = mean(true_D[si]),
    catch_reference = mean(rowSums(catch_hist[ci, c("W", "OW"), drop = FALSE])),   # 含 TR 一次 | includes TR once
    change_base_catch = sum(catch_hist[n_hist, c("W", "OW")]),
    last_total_catch = sum(catch_hist[n_hist, ]),
    initial_true_D = mean(tail(true_D, cycle)),
    catch_history = tail(catch_hist, cycle), depletion_history = tail(est_D, cycle),
    true_D = true_D, est_D = est_D)
}
cat("replicate 1: 2017–2019 estimated depletion", round(init[[1]]$state_reference, 6),
    "| 2020–2022 W+OW catch", round(init[[1]]$catch_reference, 3), "\n")

# ---- 5. 20 年闭环预测：24 个情景 × 4 个 MP × N_REP 个重复 | 20-year closed loop: 24 scenarios × 4 MPs × N_REP replicates ----
# HCR、母船队配额分配、措施选择与评估指标都与 06 相同；捕捞、观测与估计按五种渔具计算。
# HCR, parent-fleet allocation, measure choice and metrics are as in 06; fishing, observation and estimation use five gears.
# TR 是 W 母配额内的固定子配额（历史比例 tr_alloc）；没用完的 LL 配额不转给 TR。
# TR holds a fixed sub-quota inside the parent W quota (historical fraction tr_alloc); unused LL quota does not move to TR.
x1 <- hcr_points[1, 1]; y1 <- hcr_points[1, 2]           # (0.37, 0.20)
x2 <- hcr_points[2, 1]; y2 <- hcr_points[2, 2]           # (0.94, 1.09)，到 1.29 保持 1.09 | flat to 1.29
x3 <- hcr_points[3, 1]; x4 <- hcr_points[4, 1]; y4 <- hcr_points[4, 2]   # (1.29, 1.09) → (1.59, 1.29)

rows <- list(); traj1 <- list()
viol <- c(tac = 0, gear = 0, capacity = 0, cap = 0, parent = 0)   # 约束检查的最大偏差 | largest constraint deviations
for (cell_id in design$case_id) {
  cell <- design[design$case_id == cell_id, ]
  S1 <- cell$S1 == 1; S2 <- cell$S2; S3 <- cell$S3 == 1; S4 <- cell$S4 == 1
  bias_centre <- abs(misreporting) * S2                  # S2：−0.25、0 或 +0.25 | S2: −0.25, 0 or +0.25
  for (mp in mps) {
    mp_i <- match(mp, mps)                               # 1 = MP-A, 2 = MP-B, 3 = MP-C, 4 = MP-D
    for (r in seq_len(N_REP)) {
      d <- draws[[r]]; i0 <- init[[r]]
      # ---- 5.1 预测开始时的状态（同 06）| state at the start of the projection (as in 06) ----
      N <- i0$N; U <- i0$U; est <- i0$est; U_history <- i0$U_history
      target_D <- target_relative * i0$true_state_reference
      recent_catch <- i0$catch_history; recent_D <- i0$depletion_history
      nW <- round(n_segments * effort_end[["OW"]] / sum(effort_end[c("OW", "OI")]))
      seg_measure <- c(rep("W", nW), rep("I", n_segments - nW))
      seg_capacity <- rep(sum(effort_end[c("OW", "OI")]) / n_segments, n_segments)
      seg_last <- rep(-min_measure_years, n_segments)
      external_capacity <- effort_end[["E"]]
      tac <- i0$change_base_catch; first_tac <- tac
      targets <- setNames(rep(Inf, 4), fleets)
      enso <- 0; deduction <- 0
      p_WI <- 0; p_IW <- 0
      yv <- list()
      for (y in seq_len(n_proj)) {
        k <- n_hist + y; q <- q_end * effort_creep^y
        # ---- 5.2 环境状态（S3）与 E 容量增长（S4）| environmental state (S3) and E-capacity growth (S4) ----
        if (S3) enso <- c(-1, 0, 1)[which(d$enso_u[y] <= cumsum(env_transition[match(enso, -1:1), ]))[1]]
        if (S4) external_capacity <- external_capacity * (1 + epo_growth)
        capacity <- c(W = effort_end[["W"]], OW = sum(seg_capacity[seg_measure == "W"]),
                      OI = sum(seg_capacity[seg_measure == "I"]), E = external_capacity)
        cap <- if (mp_i == 4) overlap_effort_cap else Inf
        ch_WI <- 0L; ch_IW <- 0L
        decision_year <- (y - 1) %% cycle == 0
        if (decision_year) {
          # ---- 5.3 HCR 与变化限制（同 06）| HCR and change limit (as in 06) ----
          decision_D <- mean(recent_D)
          rs <- decision_D / i0$state_reference
          zz <- pmin(1, pmax(0, (rs - x1) / (x2 - x1)))
          rise <- -expm1(-hcr_curvature * (x2 - x1) * zz) / -expm1(-hcr_curvature * (x2 - x1))
          scalar <- ifelse(rs < x2, y1 + (y2 - y1) * rise,
                           y2 + (y4 - y2) * pmin(1, pmax(0, (rs - x3) / (x4 - x3))))
          hcr_tac <- i0$catch_reference * scalar
          previous_tac <- tac
          lower_lim <- change_bounds[1] * previous_tac; upper_lim <- change_bounds[2] * previous_tac
          tac <- min(upper_lim, max(lower_lim, hcr_tac))
          cl_lower <- hcr_tac < lower_lim; cl_upper <- hcr_tac > upper_lim
          # ---- 5.4 四个 MP 的母船队配额（同 06；W 配额含 TR）| parent-fleet quotas of the four MPs (as in 06; W includes TR) ----
          targets <- setNames(rep(Inf, 4), fleets); deduction <- 0
          if (mp_i == 2) {
            obs_sh <- colSums(recent_catch)
            obs_sh <- if (sum(obs_sh) > 0) obs_sh / sum(obs_sh) else share_ref
            targets <- tac / managed_share * obs_sh
          } else {
            targets[c("W", "OW")] <- tac * share_ref[c("W", "OW")] / managed_share
            if (mp_i == 3) {
              external <- mean(rowSums(recent_catch[, c("OI", "E"), drop = FALSE]))
              deduction <- max(0, external - previous_tac / managed_share * sum(share_ref[c("OI", "E")]))
              total <- sum(targets[c("W", "OW")]); deduction <- min(deduction, total)
              targets[c("W", "OW")] <- if (total > 0) targets[c("W", "OW")] * (1 - deduction / total) else 0
            }
          }
          if (mp_i == 2) {
            tot <- sum(targets[c("OW", "OI")]); den <- sum(capacity[c("OW", "OI")])
            targets[c("OW", "OI")] <- if (den > 0) tot * capacity[c("OW", "OI")] / den else 0
          }
          # ---- 5.5 预测重疊区配额利用率（估计种群、平均移动）| predicted overlap quota use (estimated numbers, mean movement) ----
          # 换措施的概率只用 OW、OI 的预测努力量，所以只解 R2 | switching uses only the OW and OI forecasts, so only R2 is solved
          Nm_p <- t(movement) %*% est
          available <- capacity; ov <- sum(available[c("OW", "OI")])
          if (ov > cap) available[c("OW", "OI")] <- available[c("OW", "OI")] * cap / ov
          Fp <- setNames(numeric(4), fleets)
          j <- 2:3; Fcap <- q * available[j]; Tj <- targets[j]; Na <- Nm_p[2, ]
          if (!(sum(Fcap) == 0 || sum(Na * w) == 0 || all(Tj == 0))) {
            req <- function(x) { Zx <- M + x * sel; pmin(Fcap, Tj / sum(sel / Zx * (-expm1(-Zx)) * Na * w)) }
            tot <- if (sum(req(sum(Fcap))) - sum(Fcap) >= -1e-12) sum(Fcap) else
              uniroot(function(x) sum(req(x)) - x, c(0, sum(Fcap)), tol = 1e-11)$root
            Fp[j] <- req(tot)
          }
          pred_effort <- Fp / q
          # ---- 5.6 S1：重疊区船队选择适用措施（同 06）| S1: overlap fleets choose applicable measures (as in 06) ----
          use_OW <- if (capacity[["OW"]] > 0) pred_effort[["OW"]] / capacity[["OW"]] else as.numeric(targets[["OW"]] > 0)
          use_OI <- if (capacity[["OI"]] > 0) pred_effort[["OI"]] / capacity[["OI"]] else as.numeric(targets[["OI"]] > 0)
          gap <- use_OI - use_OW
          p_WI <- if (gap <= switch_cost) 0 else plogis(switch_intercept + switch_slope * (gap - switch_cost))
          p_IW <- if (-gap <= switch_cost) 0 else plogis(switch_intercept + switch_slope * (-gap - switch_cost))
          if (S1) {
            eligible <- y - seg_last >= min_measure_years
            pr <- ifelse(seg_measure == "W", p_WI, p_IW)
            ch <- eligible & d$measure_change_u[y, ] < switch_strength * pr
            old <- seg_measure
            if (any(ch)) { seg_measure[ch] <- ifelse(old[ch] == "W", "I", "W"); seg_last[ch] <- as.integer(y) }
            ch_WI <- sum(ch & old == "W"); ch_IW <- sum(ch & old == "I")
            capacity <- c(W = effort_end[["W"]], OW = sum(seg_capacity[seg_measure == "W"]),
                          OI = sum(seg_capacity[seg_measure == "I"]), E = external_capacity)
            if (mp_i == 2) {
              tot <- sum(targets[c("OW", "OI")]); den <- sum(capacity[c("OW", "OI")])
              targets[c("OW", "OI")] <- if (den > 0) tot * capacity[c("OW", "OI")] / den else 0
            }
          }
        }
        # ---- 5.7 实施：W 母配额拆成 LL 与 TR 子配额，再解五种渔具的 F | implementation: W quota split into LL and TR sub-quotas, then F of five gears ----
        impl <- exp(sigma_impl * d$implementation_z[y, ] - .5 * sigma_impl^2)
        controlled <- is.finite(targets)
        attempted <- targets; attempted[controlled] <- targets[controlled] * impl[controlled]   # TR 共用 W 的实施乘子 | TR shares W's multiplier
        Py <- movement                                                    # 本年移动矩阵（S3）| this year's movement (S3)
        Py[1, 2] <- Py[1, 2] + movement_shift[1] * enso; Py[1, 1] <- Py[1, 1] - movement_shift[1] * enso
        Py[2, 3] <- Py[2, 3] + movement_shift[2] * enso; Py[2, 2] <- Py[2, 2] - movement_shift[2] * enso
        Nm <- t(Py) %*% N
        available <- capacity; ov <- sum(available[c("OW", "OI")])
        if (ov > cap) available[c("OW", "OI")] <- available[c("OW", "OI")] * cap / ov
        T5 <- c(W_LL = attempted[["W"]] * (1 - tr_alloc), attempted[c("OW", "OI", "E")], TR = attempted[["W"]] * tr_alloc)
        Fg <- setNames(numeric(5), gear_names)
        # R1 联合求解（同第 4 节）：LL 受容量 q × 可用努力量约束，TR 受 F 上限约束 | R1 joint solve (as in section 4): LL bounded by q × capacity, TR by its F ceiling
        Na <- Nm[1, ]; Tj <- T5[c("W_LL", "TR")]; Fcap <- c(q * available[["W"]], TR_Fcap)
        if (sum(Na * w) == 0 || all(Tj == 0) || all(Fcap == 0)) {
          Fg[c(1, 5)] <- Fcap * 0
        } else if (same_sel) {
          req <- function(x) { Zx <- M + x * sel; pmin(Fcap, Tj / sum(sel / Zx * (-expm1(-Zx)) * Na * w)) }
          tot <- if (sum(req(sum(Fcap))) - sum(Fcap) >= -1e-12) sum(Fcap) else
            uniroot(function(x) sum(req(x)) - x, c(0, sum(Fcap)), tol = 1e-11)$root
          Fg[c(1, 5)] <- req(tot)
        } else {
          Nw <- Na * w
          res <- function(x) {                           # 牛顿法与线搜索需要残差函数 | Newton and its line search need the residual
            Zx <- M + as.numeric(x %*% S_R1)
            gx <- as.numeric(S_R1 %*% (Nw * (-expm1(-Zx)) / Zx))
            x - pmin(Fcap, ifelse(gx > 0, Tj / gx, ifelse(Tj > 0, Inf, 0)))
          }
          g0 <- as.numeric(S_R1 %*% (Nw * (-expm1(-M)) / M))
          F2 <- pmin(Fcap, ifelse(g0 > 0, Tj / g0, Inf)); F2[Tj == 0] <- 0
          for (it in 1:60) {
            r0 <- res(F2)
            if (max(abs(r0)) < 1e-12) break
            Zx <- M + as.numeric(F2 %*% S_R1); gx <- as.numeric(S_R1 %*% (Nw * (-expm1(-Zx)) / Zx))
            free <- ifelse(gx > 0, Tj / gx, ifelse(Tj > 0, Inf, 0)) < Fcap & Tj > 0 & is.finite(Tj)
            dg <- (S_R1 * rep(Nw * ((Zx + 1) * exp(-Zx) - 1) / Zx^2, each = 2)) %*% t(S_R1)
            J <- diag(2)
            if (any(free)) J[free, ] <- J[free, , drop = FALSE] + dg[free, , drop = FALSE] * (Tj[free] / gx[free]^2)
            step <- as.numeric(solve(J, r0))
            hv <- 0
            while (max(abs(res(pmax(0, pmin(Fcap, F2 - step * 2^(-hv)))))) >= max(abs(r0))) {
              hv <- hv + 1; if (hv > 25) stop("Joint gear solver line search failed")
            }
            F2 <- pmax(0, pmin(Fcap, F2 - step * 2^(-hv)))
            if (it == 60) stop("Joint gear solver did not converge")
          }
          Fg[c(1, 5)] <- F2
        }
        for (a in 2:3) {       # R2、R3 与 06 相同 | R2 and R3 as in 06
          j <- which(gear_area == a); Fcap <- q * available[j]; Tj <- T5[j]; Na <- Nm[a, ]
          if (sum(Fcap) == 0 || sum(Na * w) == 0 || all(Tj == 0)) { Fg[j] <- Fcap * 0; next }
          req <- function(x) { Zx <- M + x * sel; pmin(Fcap, Tj / sum(sel / Zx * (-expm1(-Zx)) * Na * w)) }
          tot <- if (sum(req(sum(Fcap))) - sum(Fcap) >= -1e-12) sum(Fcap) else
            uniroot(function(x) sum(req(x)) - x, c(0, sum(Fcap)), tol = 1e-11)$root
          Fg[j] <- req(tot)
        }
        # 努力量只对 4 个延绳钓段定义；TR 由渔获控制，不用 q 换算努力量 | effort exists for the 4 LL segments only; TR is catch-controlled
        effort <- setNames(Fg[1:4] / q, fleets)
        # ---- 5.8 真实种群一年 | one year of the true population ----
        rmult <- exp(sigma_R * d$recruitment_z[k] - .5 * sigma_R^2)
        SB <- as.numeric(Nm %*% (mat * w)); VB <- as.numeric(Nm %*% (sel * w))
        Z <- rbind(Fg[1] * sel + Fg[5] * tr_sel, (Fg[2] + Fg[3]) * sel, Fg[4] * sel) + Mmat
        # 各渔具各年龄的渔获尾数 = 选择性 × F × N (1 − e^−Z) / Z | catch numbers by gear and age
        numbers <- gear_sel * (Nm * (-expm1(-Z)) / Z)[gear_area, ] * Fg
        catch_gear <- setNames(as.numeric(numbers %*% w), gear_names)
        catch <- c(W = catch_gear[["W_LL"]] + catch_gear[["TR"]], catch_gear[c("OW", "OI", "E")])
        rec <- 4 * h * R0 * sum(SB) / (SB0 * (1 - h) + sum(SB) * (5 * h - 1)) * rmult
        surv <- Nm * exp(-Z)
        N <- matrix(0, 3, A); N[, 2:A] <- surv[, 1:(A - 1)]; N[, A] <- N[, A] + surv[, A]; N[, 1] <- rec * recruit_share
        Um <- t(Py) %*% U                                                 # 未捕捞配对种群 | paired unfished population
        SBu <- as.numeric(Um %*% (mat * w))
        Zu <- outer(c(0, 0, 0), sel) + Mmat
        recu <- 4 * h * R0 * sum(SBu) / (SB0 * (1 - h) + sum(SBu) * (5 * h - 1)) * rmult
        survu <- Um * exp(-Zu)
        U <- matrix(0, 3, A); U[, 2:A] <- survu[, 1:(A - 1)]; U[, A] <- U[, A] + survu[, A]; U[, 1] <- recu * recruit_share
        # ---- 5.9 观测：S2 只改重疊区上报；TR 与 W 共用观测乘子 | observation: S2 changes overlap reports only; TR shares W's multiplier ----
        bias <- if (S2 != 0) bias_centre * (.5 + d$reporting_u[k]) else 0
        reported <- catch * exp(sigma_catch * d$catch_z[k, ] - .5 * sigma_catch^2)
        reported[c("OW", "OI")] <- reported[c("OW", "OI")] * (1 + bias)
        rep_TR <- catch_gear[["TR"]] * exp(sigma_catch * d$catch_z[k, 1] - .5 * sigma_catch^2)
        gear_rep <- c(W_LL = reported[["W"]] - rep_TR, reported[c("OW", "OI", "E")], TR = rep_TR)
        index <- VB * exp(sigma_index * d$index_z[k, ] - .5 * sigma_index^2)
        area_catch <- c(sum(reported[fleet_area == 1]), sum(reported[fleet_area == 2]), sum(reported[fleet_area == 3]))
        # ---- 5.10 估计方法：R1 用分渔具上报渔获联合反推（同第 4 节），R2、R3 同 06 | EM: R1 joint inversion of gear reports (section 4), R2/R3 as in 06 ----
        em_enso <- if (estimation_movement == "observed") enso else 0
        Pe <- movement
        Pe[1, 2] <- Pe[1, 2] + movement_shift[1] * em_enso; Pe[1, 1] <- Pe[1, 1] - movement_shift[1] * em_enso
        Pe[2, 3] <- Pe[2, 3] + movement_shift[2] * em_enso; Pe[2, 2] <- Pe[2, 2] - movement_shift[2] * em_enso
        Nm_e <- t(Pe) %*% est
        VB_pred <- pmax(as.numeric(Nm_e %*% (sel * w)), 1e-12)
        filtered <- Nm_e * exp(index_gain * log(pmax(index, 1e-12) / VB_pred))
        SB_e <- as.numeric(filtered %*% (mat * w))
        F_e <- setNames(numeric(5), gear_names)
        Na <- filtered[1, ]; Tj <- gear_rep[c("W_LL", "TR")]; Fcap <- c(em_Fmax, TR_Fcap)
        if (sum(Na * w) == 0 || all(Tj == 0) || all(Fcap == 0)) {
          F_e[c(1, 5)] <- Fcap * 0
        } else if (same_sel) {
          req <- function(x) { Zx <- M + x * sel; pmin(Fcap, Tj / sum(sel / Zx * (-expm1(-Zx)) * Na * w)) }
          tot <- if (sum(req(sum(Fcap))) - sum(Fcap) >= -1e-12) sum(Fcap) else
            uniroot(function(x) sum(req(x)) - x, c(0, sum(Fcap)), tol = 1e-11)$root
          F_e[c(1, 5)] <- req(tot)
        } else {
          Nw <- Na * w
          res <- function(x) {                           # 牛顿法与线搜索需要残差函数 | Newton and its line search need the residual
            Zx <- M + as.numeric(x %*% S_R1)
            gx <- as.numeric(S_R1 %*% (Nw * (-expm1(-Zx)) / Zx))
            x - pmin(Fcap, ifelse(gx > 0, Tj / gx, ifelse(Tj > 0, Inf, 0)))
          }
          g0 <- as.numeric(S_R1 %*% (Nw * (-expm1(-M)) / M))
          F2 <- pmin(Fcap, ifelse(g0 > 0, Tj / g0, Inf)); F2[Tj == 0] <- 0
          for (it in 1:60) {
            r0 <- res(F2)
            if (max(abs(r0)) < 1e-12) break
            Zx <- M + as.numeric(F2 %*% S_R1); gx <- as.numeric(S_R1 %*% (Nw * (-expm1(-Zx)) / Zx))
            free <- ifelse(gx > 0, Tj / gx, ifelse(Tj > 0, Inf, 0)) < Fcap & Tj > 0 & is.finite(Tj)
            dg <- (S_R1 * rep(Nw * ((Zx + 1) * exp(-Zx) - 1) / Zx^2, each = 2)) %*% t(S_R1)
            J <- diag(2)
            if (any(free)) J[free, ] <- J[free, , drop = FALSE] + dg[free, , drop = FALSE] * (Tj[free] / gx[free]^2)
            step <- as.numeric(solve(J, r0))
            hv <- 0
            while (max(abs(res(pmax(0, pmin(Fcap, F2 - step * 2^(-hv)))))) >= max(abs(r0))) {
              hv <- hv + 1; if (hv > 25) stop("Joint gear solver line search failed")
            }
            F2 <- pmax(0, pmin(Fcap, F2 - step * 2^(-hv)))
            if (it == 60) stop("Joint gear solver did not converge")
          }
          F_e[c(1, 5)] <- F2
        }
        for (a in 2:3) {
          jj <- c(2, 4)[a - 1]
          Na <- filtered[a, ]; Ta <- area_catch[a]
          if (sum(Na * w) == 0 || Ta == 0) next
          req <- function(x) { Zx <- M + x * sel; pmin(em_Fmax, Ta / sum(sel / Zx * (-expm1(-Zx)) * Na * w)) }
          tot <- if (sum(req(em_Fmax)) - em_Fmax >= -1e-12) em_Fmax else
            uniroot(function(x) sum(req(x)) - x, c(0, em_Fmax), tol = 1e-11)$root
          F_e[jj] <- req(tot)
        }
        Ze <- rbind(F_e[1] * sel + F_e[5] * tr_sel, (F_e[2] + F_e[3]) * sel, F_e[4] * sel) + Mmat
        rec_e <- 4 * h * R0 * sum(SB_e) / (SB0 * (1 - h) + sum(SB_e) * (5 * h - 1))
        surve <- filtered * exp(-Ze)
        est <- matrix(0, 3, A); est[, 2:A] <- surve[, 1:(A - 1)]; est[, A] <- est[, A] + surve[, A]; est[, 1] <- rec_e * recruit_share
        D_est <- sum(SB_e) / SB0
        recent_catch <- rbind(recent_catch[-1, , drop = FALSE], reported)   # 母船队上报渔获 | parent-fleet reports
        recent_D <- c(tail(recent_D, cycle - 1), D_est)
        # ---- 5.11 本年记录（06 的量 + 渔具量）| this year's record (quantities of 06 + gear quantities) ----
        D <- sum(SB) / mean(tail(U_history, unfished_window))
        U_history <- c(tail(U_history, unfished_window - 1), sum(SBu))
        C <- sum(catch); finite_sum <- sum(targets[controlled])
        uncapped <- !controlled & !c(FALSE, is.finite(cap), is.finite(cap), FALSE)
        finite_gear <- is.finite(T5)
        yv[[y]] <- list(
          year = y, D = D, D_true_fixed = sum(SB) / SB0, D_est = D_est, D_overlap = SB[2] / SBu[2],
          decision_year = decision_year, cl_lower = cl_lower, cl_upper = cl_upper, rs = rs,
          hcr_tac = hcr_tac, tac = tac, C_total = C,
          not_managed = if (C > 0) sum(catch[uncapped]) / C else 0,
          outside_W = if (C > 0) sum(catch[c("OI", "E")]) / C else 0,
          exceed = if (finite_sum > 0) sum(pmax(catch[controlled] - targets[controlled], 0)) / finite_sum else 0,
          uncaught = if (finite_sum > 0) sum(pmax(targets[controlled] - catch[controlled], 0)) / finite_sum else 0,
          tac_violation = max(c(0, catch[controlled] - attempted[controlled])),
          capacity_violation = max(c(0, effort - available)),
          E_overlap = sum(effort[c("OW", "OI")]),
          cap_violation = if (is.finite(cap)) max(0, sum(effort[c("OW", "OI")]) - cap) else 0,
          changes = ch_WI + ch_IW, deduction = deduction,
          catch = catch, attempted = ifelse(controlled, attempted, NA_real_),
          # 渔具量 | gear quantities
          C_TR = catch_gear[["TR"]], C_W_LL = catch_gear[["W_LL"]], C_LL = sum(catch_gear[1:4]),
          F_TR = Fg[["TR"]], T_TR = T5[["TR"]], T_WLL = T5[["W_LL"]],
          TR_bound = Fg[["TR"]] >= TR_Fmax * (1 - 1e-8),                 # TR 是否触及 F 上限 | TR at its F ceiling
          gear_violation = max(c(0, catch_gear[finite_gear] - T5[finite_gear])))
      }
      # ---- 5.12 汇总一个重复：06 的指标 + 渔具指标（评估期 = 第 11–20 年）| summarise one replicate: metrics of 06 + gear metrics ----
      Dv <- sapply(yv, `[[`, "D"); Cv <- sapply(yv, `[[`, "C_total")
      Cf <- t(sapply(yv, `[[`, "catch")); Af <- t(sapply(yv, `[[`, "attempted"))
      e <- evaluation_years
      De <- Dv[e]; Ce <- Cv[e]
      Dt <- sapply(yv, `[[`, "D_true_fixed")[e]; Dm <- sapply(yv, `[[`, "D_est")[e]
      dec <- sapply(yv, `[[`, "decision_year")
      if (i0$initial_true_D >= target_D) recovered <- c(1, 0) else {   # 3 年滑动平均达到目标 | 3-year running mean reaches target
        m3 <- numeric(0)
        for (t3 in 3:n_proj) m3 <- c(m3, mean(Dv[(t3 - 2):t3]))
        hit <- which(m3 >= target_D)
        recovered <- if (length(hit)) c(1, hit[1] + 2) else c(0, n_proj)
      }
      full_caught <- rep(NA_real_, 4)                                # 配额被捕满的年份比例 | share of years quota fully caught
      for (fi in 1:4) {
        qa <- Af[e, fi]; act <- is.finite(qa) & qa > 0
        if (any(act)) full_caught[fi] <- mean(Cf[e, fi][act] >= qa[act] * (1 - 1e-6))
      }
      # 渔具指标 | gear metrics
      CTR <- sapply(yv, `[[`, "C_TR")[e]; CWLL <- sapply(yv, `[[`, "C_W_LL")[e]; CLL <- sapply(yv, `[[`, "C_LL")[e]
      TTR <- sapply(yv, `[[`, "T_TR")[e]; TWLL <- sapply(yv, `[[`, "T_WLL")[e]
      act_TR <- is.finite(TTR) & TTR > 0; act_WLL <- is.finite(TWLL) & TWLL > 0   # 只计配额为正的年份 | positive-quota years only
      rows[[length(rows) + 1]] <- data.frame(
        case_id = cell_id, mp = mp, rep = r,
        P_below_LRP = mean(De < limit_reference), P_ever_LRP = as.numeric(any(De < limit_reference)),
        P_below_42 = mean(De < diagnostic_threshold), P_ever_42 = as.numeric(any(De < diagnostic_threshold)),
        mean_D = mean(De), mean_relative_D = mean(De / i0$true_state_reference),
        mae_target = mean(abs(De - target_D)),
        mean_target = mean(rep(target_D, length(e))), signed_target_gap = mean(De - target_D),
        mean_catch = mean(Ce),
        AAV = if (sum(head(Ce, -1)) > 0) sum(abs(diff(Ce))) / sum(head(Ce, -1)) else NA_real_,
        not_managed_fraction = mean(sapply(yv, `[[`, "not_managed")[e]),
        outside_W_fraction = mean(sapply(yv, `[[`, "outside_W")[e]),
        D_overlap = mean(sapply(yv, `[[`, "D_overlap")[e]),
        overlap_effort = mean(sapply(yv, `[[`, "E_overlap")[e]),
        quota_exceedance_ratio = mean(sapply(yv, `[[`, "exceed")[e]),
        uncaught_quota_ratio = mean(sapply(yv, `[[`, "uncaught")[e]),
        effort_cap_exceedance = mean(sapply(yv, `[[`, "cap_violation")[e] > 1e-8),
        measure_changes_in_window = sum(sapply(yv, `[[`, "changes")[e]),
        mean_deduction = mean(sapply(yv, `[[`, "deduction")[e]),
        catch_W = mean(Cf[e, 1]), catch_OW = mean(Cf[e, 2]), catch_OI = mean(Cf[e, 3]), catch_E = mean(Cf[e, 4]),
        change_limit_lower_frequency = mean(sapply(yv, `[[`, "cl_lower")[dec]),
        change_limit_upper_frequency = mean(sapply(yv, `[[`, "cl_upper")[dec]),
        change_base_initial_cut = max(0, 1 - first_tac / i0$change_base_catch),
        first_tac_cut = max(0, 1 - yv[[1]]$tac / i0$change_base_catch),
        early_catch_reduction = 1 - mean(Cv[1:3]) / i0$last_total_catch,
        recovery_event = recovered[1], recovery_time_restricted = recovered[2],
        relative_catch = mean(Ce) / i0$last_total_catch,
        P_below_target = mean(De < target_D),
        mean_D_true_fixed = mean(Dt), mean_D_est_fixed = mean(Dm),
        estimation_bias_fixed = mean(Dm - Dt), estimation_rmse_fixed = sqrt(mean((Dm - Dt)^2)),
        mean_denominator_gap = mean(De - Dt),
        quota_fully_caught_W = full_caught[1], quota_fully_caught_OW = full_caught[2],
        quota_fully_caught_OI = full_caught[3], quota_fully_caught_E = full_caught[4],
        mean_TR_catch = mean(CTR), mean_LL_catch = mean(CLL),
        TR_weight_fraction = sum(CTR) / sum(Ce),                      # TR 占总渔获重量 | TR share of total catch weight
        quota_fully_caught_TR = if (!any(act_TR)) NA_real_ else mean(CTR[act_TR] >= TTR[act_TR] * (1 - 1e-6)),
        quota_fully_caught_W_LL = if (!any(act_WLL)) NA_real_ else mean(CWLL[act_WLL] >= TWLL[act_WLL] * (1 - 1e-6)),
        TR_F_bound_frequency = mean(sapply(yv, `[[`, "TR_bound")[e]),
        S1 = cell$S1, S2 = cell$S2, S3 = cell$S3, S4 = cell$S4, config_id = GEAR_CONFIG)
      viol <- pmax(viol, c(max(sapply(yv, `[[`, "tac_violation")), max(sapply(yv, `[[`, "gear_violation")),
                           max(sapply(yv, `[[`, "capacity_violation")), max(sapply(yv, `[[`, "cap_violation")),
                           max(abs(Cv - sapply(yv, `[[`, "C_LL") - sapply(yv, `[[`, "C_TR")))))
      # 保存第 1 个重复的逐年结果（用于第 6 节）| keep replicate 1 year by year (for section 6)
      if (r == 1) traj1[[length(traj1) + 1]] <- data.frame(
        case_id = cell_id, mp = mp, year = seq_len(n_proj), decision = dec,
        rel_status = round(sapply(yv, `[[`, "rs"), 3), TAC = round(sapply(yv, `[[`, "tac"), 1),
        W_LL = round(sapply(yv, `[[`, "C_W_LL"), 1), TR = round(sapply(yv, `[[`, "C_TR"), 1),
        OW = round(Cf[, 2], 1), OI = round(Cf[, 3], 1), E = round(Cf[, 4], 1),
        F_TR = round(sapply(yv, `[[`, "F_TR"), 4), depletion = round(Dv, 4))
    }
  }
  cat(format(Sys.time(), "%H:%M:%S"), cell_id, "done\n")
}
metrics <- bind_rows(rows)
# 数值检查：母配额、子配额、容量与上限都没有被突破，W = LL + TR | numerical checks: quotas, capacity and cap hold; W = LL + TR
cat("largest constraint deviations:\n"); print(signif(viol, 3))
stopifnot(viol[["tac"]] < 1e-7, viol[["gear"]] < 1e-7, viol[["capacity"]] < 1e-10, viol[["cap"]] < 1e-10,
          viol[["parent"]] < 1e-7, all(is.finite(metrics$mean_D)))

# ---- 6. 看结果 | look at the results ----
# 6a. 情景 C05（S1–S4 全关）第 1 个重复的逐年结果 | replicate 1 of scenario C05 (S1–S4 off), year by year
tr1 <- bind_rows(traj1)
for (mp in mps) {
  cat("\n====", mp, "| scenario C05 | replicate 1 ====\n")
  print(tr1[tr1$case_id == "C05" & tr1$mp == mp, -(1:2)], row.names = FALSE)
}
# 6b. C05 中 TR 的表现 | TR outcomes in C05
cat("\nTR outcomes in scenario C05 (mean over replicates):\n")
print(metrics |> filter(case_id == "C05") |> group_by(mp) |>
        summarise(TR_catch = mean(mean_TR_catch), TR_weight_fraction = mean(TR_weight_fraction),
                  TR_quota_fully_caught = mean(quota_fully_caught_TR), TR_F_at_bound = mean(TR_F_bound_frequency),
                  .groups = "drop"), width = Inf)
# 6c. 综合情景 C23（S1 开、S2 = +1、S3 开）相对 C05 的配对效应 | paired effect of the combined scenario C23 (S1 on, S2 = +1, S3 on) vs C05
paired <- inner_join(metrics |> filter(case_id == "C05") |> select(mp, rep, D_C05 = mean_D, catch_C05 = relative_catch),
                     metrics |> filter(case_id == "C23") |> select(mp, rep, D_C23 = mean_D, catch_C23 = relative_catch),
                     by = c("mp", "rep"))
cat("\ncombined scenario C23 vs C05 (paired within replicate, then averaged):\n")
print(paired |> group_by(mp) |>
        summarise(n = n(), depletion_C05 = mean(D_C05), depletion_C23 = mean(D_C23), effect_depletion = mean(D_C23 - D_C05),
                  rel_catch_C05 = mean(catch_C05), rel_catch_C23 = mean(catch_C23), effect_rel_catch = mean(catch_C23 - catch_C05),
                  .groups = "drop"), width = Inf)
cat("elapsed:", format(Sys.time() - t0), "\n")

# ---- 7. 保存 | save ----
con <- gzfile(file.path(dest, "replicate_metrics.csv.gz"), "wt")
write.csv(metrics, con, row.names = FALSE, na = ""); close(con)
cat("10 done:", dest, "\n")
