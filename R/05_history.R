# 05_history.R — 随机数与 30 年合成历史（每个重复一条）：真实种群、观测、估计
# 05_history.R — random numbers and the 30-year synthetic history (one per replicate): true population,
#                observations, estimation
#
# 历史阶段结束时得到 HCR 需要的三个参照量 | the history ends with the three references the HCR needs:
#   state_reference  2017–2019 估计耗竭平均（相对状态的分母）| mean estimated depletion (denominator of relative status)
#   catch_reference  2020–2022 W+OW 观测渔获平均（渔获乘数的基数）| mean observed W+OW catch (base of the catch scalar)
#   change_base_catch 2024 年 W+OW 观测渔获（第一次 −5%/+10% 的比较基数）| 2024 W+OW catch (base of the first change limit)

if (!exists("q_end")) source("R/04_initialise.R")

# ---- 1. 预先抽取全部随机数（每个重复一块，顺序固定）| draw all random numbers in advance (one block per replicate, fixed order) ----
# 重复 r 的随机数与情景、MP 无关，因此同一重复可在情景和 MP 之间配对比较。
# Replicate r's numbers do not depend on scenario or MP, so replicates pair across scenarios and MPs.
set.seed(SEED)
n_all <- n_hist + n_proj
draws <- vector("list", N_REP)
for (r in seq_len(N_REP)) {
  draws[[r]] <- list(
    recruitment_z    = rnorm(n_all),                                   # 补充偏差 | recruitment deviations
    catch_z          = matrix(rnorm(n_all * 4), n_all, 4),             # 渔获观测误差 | catch observation errors
    index_z          = matrix(rnorm(n_all * 3), n_all, 3),             # 指数误差 | index errors
    implementation_z = matrix(rnorm(n_proj * 4), n_proj, 4),           # 实施误差 | implementation errors
    enso_u           = runif(n_proj),                                  # S3 环境状态 | S3 environmental state
    measure_change_u = matrix(runif(n_proj * n_segments), n_proj, n_segments),  # S1 更换措施 | S1 changes
    reporting_u      = runif(n_all))                                   # S2 误报幅度 | S2 misreporting size
}

# ---- 2. 每个重复跑 30 年历史 | run the 30-year history for each replicate ----
init <- vector("list", N_REP)
for (r in seq_len(N_REP)) {
  d <- draws[[r]]
  N <- N_start; est <- N_start; U <- N0                 # 真实、估计、未捕捞配对种群 | true, estimated, paired unfished
  U_history <- rep(SB0, unfished_window)                # 未捕捞产卵潜力的记录 | record of unfished spawning potential
  catch_hist <- matrix(0, n_hist, 4, dimnames = list(NULL, fleets))
  est_D <- true_D <- numeric(n_hist)
  for (yr in seq_len(n_hist)) {
    rmult <- exp(sigma_R * d$recruitment_z[yr] - .5 * sigma_R^2)       # 补充乘子（参考配置 = 1）| recruitment multiplier
    # -- 真实种群 | true population --
    q <- q_initial * effort_creep^(yr - 1)
    Fy <- q * history_effort[yr, ]
    Nm <- t(movement) %*% N
    SB <- as.numeric(Nm %*% (mat * w)); VB <- as.numeric(Nm %*% (sel * w))
    F_area <- c(sum(Fy[fleet_area == 1]), sum(Fy[fleet_area == 2]), sum(Fy[fleet_area == 3]))
    Z <- outer(F_area, sel) + Mmat
    g <- c(sum(sel / Z[1, ] * (-expm1(-Z[1, ])) * Nm[1, ] * w),
           sum(sel / Z[2, ] * (-expm1(-Z[2, ])) * Nm[2, ] * w),
           sum(sel / Z[3, ] * (-expm1(-Z[3, ])) * Nm[3, ] * w))
    catch <- setNames(Fy * g[fleet_area], fleets)
    rec <- 4 * h * R0 * sum(SB) / (SB0 * (1 - h) + sum(SB) * (5 * h - 1)) * rmult
    surv <- Nm * exp(-Z)
    N <- matrix(0, 3, A); N[, 2:A] <- surv[, 1:(A - 1)]; N[, A] <- N[, A] + surv[, A]; N[, 1] <- rec * recruit_share
    # -- 未捕捞配对种群（同一补充偏差，用作耗竭的分母）| paired unfished population (same deviations; depletion denominator) --
    Um <- t(movement) %*% U
    SBu <- as.numeric(Um %*% (mat * w))
    Zu <- outer(c(0, 0, 0), sel) + Mmat
    recu <- 4 * h * R0 * sum(SBu) / (SB0 * (1 - h) + sum(SBu) * (5 * h - 1)) * rmult
    survu <- Um * exp(-Zu)
    U <- matrix(0, 3, A); U[, 2:A] <- survu[, 1:(A - 1)]; U[, A] <- U[, A] + survu[, A]; U[, 1] <- recu * recruit_share
    # -- 观测：渔获与按可捕生物量的指数 | observations: catch and a biomass index --
    obs_catch <- catch * exp(sigma_catch * d$catch_z[yr, ] - .5 * sigma_catch^2)
    index <- VB * exp(sigma_index * d$index_z[yr, ] - .5 * sigma_index^2)
    area_catch <- c(sum(obs_catch[fleet_area == 1]), sum(obs_catch[fleet_area == 2]), sum(obs_catch[fleet_area == 3]))
    # -- 估计方法：用指数修正上年估计，再按报告渔获反推 F 并推进一年 --
    # -- estimation method: correct last year's estimate with the index, back-calculate F from reported catch, project one year --
    Nm_e <- t(movement) %*% est
    VB_pred <- pmax(as.numeric(Nm_e %*% (sel * w)), 1e-12)
    filtered <- Nm_e * exp(index_gain * log(pmax(index, 1e-12) / VB_pred))   # 逐区缩放 | scale each region
    SB_e <- as.numeric(filtered %*% (mat * w))
    F_e <- numeric(3)
    for (a in 1:3) {                                   # 每区解 Baranov：报告渔获 = F × g(F) | solve Baranov per region
      Na <- filtered[a, ]; Ta <- area_catch[a]
      if (sum(Na * w) == 0 || Ta == 0) next
      req <- function(x) { Zx <- M + x * sel; pmin(em_Fmax, Ta / sum(sel / Zx * (-expm1(-Zx)) * Na * w)) }
      tot <- if (sum(req(em_Fmax)) - em_Fmax >= -1e-12) em_Fmax else
        uniroot(function(x) sum(req(x)) - x, c(0, em_Fmax), tol = 1e-11)$root
      F_e[a] <- req(tot)
    }
    Ze <- outer(F_e, sel) + Mmat
    rec_e <- 4 * h * R0 * sum(SB_e) / (SB0 * (1 - h) + sum(SB_e) * (5 * h - 1))
    surve <- filtered * exp(-Ze)
    est <- matrix(0, 3, A); est[, 2:A] <- surve[, 1:(A - 1)]; est[, A] <- est[, A] + surve[, A]; est[, 1] <- rec_e * recruit_share
    # -- 记录 | record --
    true_D[yr] <- sum(SB) / mean(tail(U_history, unfished_window))   # 真实耗竭 | true depletion
    U_history <- c(U_history, sum(SBu))
    est_D[yr] <- sum(SB_e) / SB0                                       # 估计耗竭（固定 SB0）| estimated depletion (fixed SB0)
    catch_hist[yr, ] <- obs_catch
  }
  si <- match(state_reference_years, history_years); ci <- match(catch_reference_years, history_years)
  init[[r]] <- list(
    N = N, U = U, est = est, U_history = tail(U_history, unfished_window),
    state_reference = mean(est_D[si]),                                   # HCR 状态分母 | HCR status denominator
    true_state_reference = mean(true_D[si]),                             # 目标 = 0.96 × 此值 | target = 0.96 × this
    catch_reference = mean(rowSums(catch_hist[ci, c("W", "OW"), drop = FALSE])),
    change_base_catch = sum(catch_hist[n_hist, c("W", "OW")]),
    last_total_catch = sum(catch_hist[n_hist, ]),                        # 渔获的相对尺度 | scale of relative catch
    initial_true_D = mean(tail(true_D, cycle)),
    catch_history = tail(catch_hist, cycle), depletion_history = tail(est_D, cycle),
    true_D = true_D, est_D = est_D, catch_all = catch_hist)            # 30 年观测渔获（11 用）| all 30 years (used in step 11)
}

# ---- 3. 看第 1 个重复 | look at replicate 1 ----
i1 <- init[[1]]
cat("replicate 1: 2017–2019 estimated depletion", round(i1$state_reference, 4),
    "| target", round(target_relative * i1$true_state_reference, 4),
    "| 2020–2022 W+OW catch", round(i1$catch_reference, 2), "\n")
png(file.path(OUT, "05_history_rep1.png"), width = 1400, height = 900, res = 200)
plot(history_years, i1$true_D, type = "l", lwd = 2, ylim = c(.40, .50),
     xlab = "Year (synthetic history)", ylab = "Depletion",
     main = "Replicate 1: true (line) and estimated (points) depletion")
points(history_years, i1$est_D, pch = 16, cex = .6)
invisible(dev.off())

cat("05 done:", N_REP, "replicates\n")
