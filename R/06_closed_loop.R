# 06_closed_loop.R — 20 年闭环预测：HCR → 配额分配 → 实施 → 种群 → 观测 → 估计
# 06_closed_loop.R — 20-year closed-loop projection: HCR → allocation → implementation → population
#                    → observation → estimation
#
# 单独运行时：参考情景 C05（S1–S4 全关），四个候选 MP，逐年打印第 1 个重复。
# Run on its own: reference scenario C05 (S1–S4 off), four candidate MPs, replicate 1 printed year by year.
# 07 会把 RUN_CELLS 设为全部 24 个情景再 source 本文件。
# Step 07 sets RUN_CELLS to all 24 scenarios and sources this file.
#
# 管理者只看到观测与估计，从不读取真实种群 | the manager sees observations and estimates only, never the true population

if (!exists("init")) source("R/05_history.R")

# ---- 0. 本次运行的选择（可在 source 之前覆盖）| choices for this run (override before sourcing) ----
if (!exists("RUN_CELLS"))    RUN_CELLS <- "C05"          # 情景 | scenarios
if (!exists("RUN_MPS"))      RUN_MPS <- mps              # 候选 MP | candidate MPs
if (!exists("FIXED_SHARES")) FIXED_SHARES <- FALSE       # TRUE = MP-B 固定上报份额（补充实验，见 09）| fixed shares under MP-B (step 09)
VERBOSE <- length(RUN_CELLS) == 1 && N_REP <= 5          # 小运行才逐年打印 | print years only for small runs

# HCR 控制点 | HCR control points
x1 <- hcr_points[1, 1]; y1 <- hcr_points[1, 2]           # (0.37, 0.20)
x2 <- hcr_points[2, 1]; y2 <- hcr_points[2, 2]           # (0.94, 1.09)，到 1.29 保持 1.09 | flat to 1.29
x3 <- hcr_points[3, 1]; x4 <- hcr_points[4, 1]; y4 <- hcr_points[4, 2]   # (1.29, 1.09) → (1.59, 1.29)

rows <- list(); traj <- list(); diag_rows <- list()
for (cell_id in RUN_CELLS) {
  cell <- design[design$case_id == cell_id, ]
  S1 <- cell$S1 == 1; S2 <- cell$S2; S3 <- cell$S3 == 1; S4 <- cell$S4 == 1
  bias_centre <- abs(misreporting) * S2                  # S2：−0.25、0 或 +0.25 | S2: −0.25, 0 or +0.25
  for (mp in RUN_MPS) {
    mp_i <- match(mp, mps)                               # 1 = MP-A, 2 = MP-B, 3 = MP-C, 4 = MP-D
    diag_max <- c(tac = 0, capacity = 0, cap = 0, cap_n = 0, shortfall = 0)
    for (r in seq_len(N_REP)) {
      d <- draws[[r]]; i0 <- init[[r]]
      # ---- 1. 预测开始时的状态 | state at the start of the projection ----
      N <- i0$N; U <- i0$U; est <- i0$est; U_history <- i0$U_history
      target_D <- target_relative * i0$true_state_reference          # 评估目标 | evaluation target
      recent_catch <- i0$catch_history; recent_D <- i0$depletion_history   # 最近 3 年 | last 3 years
      # 重疊区 60 个决策单元：40 个适用 WCPFC 措施，20 个适用 IATTC 措施（容量 0.2 : 0.1）
      # 60 decision units in the overlap: 40 under WCPFC and 20 under IATTC measures (capacity 0.2 : 0.1)
      nW <- round(n_segments * effort_end[["OW"]] / sum(effort_end[c("OW", "OI")]))
      seg_measure <- c(rep("W", nW), rep("I", n_segments - nW))
      seg_capacity <- rep(sum(effort_end[c("OW", "OI")]) / n_segments, n_segments)
      seg_last <- rep(-min_measure_years, n_segments)               # 上次更换的年份 | year of last change
      external_capacity <- effort_end[["E"]]
      tac <- i0$change_base_catch; first_tac <- tac                  # 第一次变化的比较基数 | base of the first change
      targets <- setNames(rep(Inf, 4), fleets)                       # Inf = 不受配额约束 | Inf = no quota
      enso <- 0; n_changes <- 0L; deduction <- 0
      p_WI <- 0; p_IW <- 0
      if (FIXED_SHARES) fixed_shares <- colSums(i0$catch_history) / sum(i0$catch_history)
      yv <- list()                                                   # 逐年记录 | annual record
      for (y in seq_len(n_proj)) {
        k <- n_hist + y; q <- q_end * effort_creep^y
        # ---- 2. 环境状态（S3）与 E 容量增长（S4）| environmental state (S3) and E-capacity growth (S4) ----
        if (S3) enso <- c(-1, 0, 1)[which(d$enso_u[y] <= cumsum(env_transition[match(enso, -1:1), ]))[1]]
        if (S4) external_capacity <- external_capacity * (1 + epo_growth)
        capacity <- c(W = effort_end[["W"]], OW = sum(seg_capacity[seg_measure == "W"]),
                      OI = sum(seg_capacity[seg_measure == "I"]), E = external_capacity)
        cap <- if (mp_i == 4) overlap_effort_cap else Inf             # MP-D 重疊区努力量上限 | MP-D cap
        ch_WI <- 0L; ch_IW <- 0L
        decision_year <- (y - 1) %% cycle == 0                       # 第 1、4、7…19 年决策 | decisions in years 1, 4, …, 19
        if (decision_year) {
          # ---- 3. HCR：相对状态 → 渔获乘数 → 建议渔获 | HCR: relative status → catch scalar → HCR catch ----
          decision_D <- mean(recent_D)                               # 最近 3 年估计耗竭 | last 3 estimated depletions
          rs <- decision_D / i0$state_reference                      # 相对 2017–2019 | relative to 2017–2019
          zz <- pmin(1, pmax(0, (rs - x1) / (x2 - x1)))
          rise <- -expm1(-hcr_curvature * (x2 - x1) * zz) / -expm1(-hcr_curvature * (x2 - x1))
          scalar <- ifelse(rs < x2, y1 + (y2 - y1) * rise,
                           y2 + (y4 - y2) * pmin(1, pmax(0, (rs - x3) / (x4 - x3))))
          hcr_tac <- i0$catch_reference * scalar
          # ---- 4. 变化限制：上次 TAC 的 95%–110% | change limit: 95–110% of the previous TAC ----
          previous_tac <- tac
          lower_lim <- change_bounds[1] * previous_tac; upper_lim <- change_bounds[2] * previous_tac
          tac <- min(upper_lim, max(lower_lim, hcr_tac))
          cl_lower <- hcr_tac < lower_lim; cl_upper <- hcr_tac > upper_lim
          # ---- 5. 四个 MP 的配额分配 | allocation under the four MPs ----
          targets <- setNames(rep(Inf, 4), fleets); deduction <- 0
          if (mp_i == 2) {
            # MP-B：共同 TAC 覆盖四支船队，按最近 3 年上报渔获份额分配 | joint TAC over all fleets, by recent reported shares
            if (FIXED_SHARES) {
              targets <- tac / sum(share_ref[c("W", "OW")]) * fixed_shares
            } else {
              obs_sh <- colSums(recent_catch)
              obs_sh <- if (sum(obs_sh) > 0) obs_sh / sum(obs_sh) else share_ref
              targets <- tac / managed_share * obs_sh
            }
          } else {
            # MP-A、C、D：TAC 只给 W 和 OW，按固定参考份额 | TAC for W and OW only, by fixed reference shares
            targets[c("W", "OW")] <- tac * share_ref[c("W", "OW")] / managed_share
            if (mp_i == 3) {
              # MP-C：OI+E 上报渔获超出假设水平的部分，从 W、OW 配额中扣除 | deduct OI+E catch above the assumed level
              external <- mean(rowSums(recent_catch[, c("OI", "E"), drop = FALSE]))
              deduction <- max(0, external - previous_tac / managed_share * sum(share_ref[c("OI", "E")]))
              total <- sum(targets[c("W", "OW")]); deduction <- min(deduction, total)
              targets[c("W", "OW")] <- if (total > 0) targets[c("W", "OW")] * (1 - deduction / total) else 0
            }
          }
          if (mp_i == 2) {   # MP-B：重疊区配额按 OW、OI 当前容量分 | MP-B: overlap quota split by current OW/OI capacity
            tot <- sum(targets[c("OW", "OI")]); den <- sum(capacity[c("OW", "OI")])
            targets[c("OW", "OI")] <- if (den > 0) tot * capacity[c("OW", "OI")] / den else 0
          }
          # ---- 6. 预测配额利用率（用估计种群、平均移动）| predicted quota use (estimated numbers, mean movement) ----
          Nm_p <- t(movement) %*% est
          available <- capacity; ov <- sum(available[c("OW", "OI")])
          if (ov > cap) available[c("OW", "OI")] <- available[c("OW", "OI")] * cap / ov
          Fp <- setNames(numeric(4), fleets)
          for (a in 1:3) {     # 解 Baranov：配额与容量同时约束 | solve Baranov under quota and capacity
            j <- which(fleet_area == a); Fcap <- q * available[j]; Tj <- targets[j]; Na <- Nm_p[a, ]
            if (sum(Fcap) == 0 || sum(Na * w) == 0 || all(Tj == 0)) { Fp[j] <- Fcap * 0; next }
            req <- function(x) { Zx <- M + x * sel; pmin(Fcap, Tj / sum(sel / Zx * (-expm1(-Zx)) * Na * w)) }
            tot <- if (sum(req(sum(Fcap))) - sum(Fcap) >= -1e-12) sum(Fcap) else
              uniroot(function(x) sum(req(x)) - x, c(0, sum(Fcap)), tol = 1e-11)$root
            Fp[j] <- req(tot)
          }
          pred_effort <- Fp / q
          # ---- 7. S1：重疊区船队选择适用措施 | S1: overlap fleets choose applicable measures ----
          # 利用率差 = OI 利用率 − OW 利用率；差超过 0.05 时按 logit 概率换到另一方 | switch by a logit of the utilisation gap
          use_OW <- if (capacity[["OW"]] > 0) pred_effort[["OW"]] / capacity[["OW"]] else as.numeric(targets[["OW"]] > 0)
          use_OI <- if (capacity[["OI"]] > 0) pred_effort[["OI"]] / capacity[["OI"]] else as.numeric(targets[["OI"]] > 0)
          gap <- use_OI - use_OW
          p_WI <- if (gap <= switch_cost) 0 else plogis(switch_intercept + switch_slope * (gap - switch_cost))
          p_IW <- if (-gap <= switch_cost) 0 else plogis(switch_intercept + switch_slope * (-gap - switch_cost))
          if (S1) {
            eligible <- y - seg_last >= min_measure_years                 # 至少 3 年 | at least 3 years
            pr <- ifelse(seg_measure == "W", p_WI, p_IW)
            ch <- eligible & d$measure_change_u[y, ] < switch_strength * pr
            old <- seg_measure
            if (any(ch)) { seg_measure[ch] <- ifelse(old[ch] == "W", "I", "W"); seg_last[ch] <- as.integer(y) }
            ch_WI <- sum(ch & old == "W"); ch_IW <- sum(ch & old == "I")
            n_changes <- n_changes + ch_WI + ch_IW
            capacity <- c(W = effort_end[["W"]], OW = sum(seg_capacity[seg_measure == "W"]),
                          OI = sum(seg_capacity[seg_measure == "I"]), E = external_capacity)
            if (mp_i == 2) {
              tot <- sum(targets[c("OW", "OI")]); den <- sum(capacity[c("OW", "OI")])
              targets[c("OW", "OI")] <- if (den > 0) tot * capacity[c("OW", "OI")] / den else 0
            }
          }
        }
        # ---- 8. 实施：实际捕捞受配额与容量约束 | implementation: catch limited by quota and capacity ----
        impl <- exp(sigma_impl * d$implementation_z[y, ] - .5 * sigma_impl^2)
        controlled <- is.finite(targets)
        attempted <- targets; attempted[controlled] <- targets[controlled] * impl[controlled]
        Py <- movement                                                    # 本年移动矩阵（S3）| this year's movement (S3)
        Py[1, 2] <- Py[1, 2] + movement_shift[1] * enso; Py[1, 1] <- Py[1, 1] - movement_shift[1] * enso
        Py[2, 3] <- Py[2, 3] + movement_shift[2] * enso; Py[2, 2] <- Py[2, 2] - movement_shift[2] * enso
        Nm <- t(Py) %*% N
        available <- capacity; ov <- sum(available[c("OW", "OI")])
        if (ov > cap) available[c("OW", "OI")] <- available[c("OW", "OI")] * cap / ov
        Fy <- setNames(numeric(4), fleets)
        for (a in 1:3) {       # 同上的 Baranov 求解 | same Baranov solve as above
          j <- which(fleet_area == a); Fcap <- q * available[j]; Tj <- attempted[j]; Na <- Nm[a, ]
          if (sum(Fcap) == 0 || sum(Na * w) == 0 || all(Tj == 0)) { Fy[j] <- Fcap * 0; next }
          req <- function(x) { Zx <- M + x * sel; pmin(Fcap, Tj / sum(sel / Zx * (-expm1(-Zx)) * Na * w)) }
          tot <- if (sum(req(sum(Fcap))) - sum(Fcap) >= -1e-12) sum(Fcap) else
            uniroot(function(x) sum(req(x)) - x, c(0, sum(Fcap)), tol = 1e-11)$root
          Fy[j] <- req(tot)
        }
        effort <- Fy / q
        # ---- 9. 真实种群一年 | one year of the true population ----
        rmult <- exp(sigma_R * d$recruitment_z[k] - .5 * sigma_R^2)
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
        # 未捕捞配对种群 | paired unfished population
        Um <- t(Py) %*% U
        SBu <- as.numeric(Um %*% (mat * w))
        Zu <- outer(c(0, 0, 0), sel) + Mmat
        recu <- 4 * h * R0 * sum(SBu) / (SB0 * (1 - h) + sum(SBu) * (5 * h - 1)) * rmult
        survu <- Um * exp(-Zu)
        U <- matrix(0, 3, A); U[, 2:A] <- survu[, 1:(A - 1)]; U[, A] <- U[, A] + survu[, A]; U[, 1] <- recu * recruit_share
        # ---- 10. 观测：S2 把重疊区上报渔获乘 (1 + b) | observation: S2 multiplies reported overlap catch by (1 + b) ----
        bias <- if (S2 != 0) bias_centre * (.5 + d$reporting_u[k]) else 0
        reported <- catch * exp(sigma_catch * d$catch_z[k, ] - .5 * sigma_catch^2)
        reported[c("OW", "OI")] <- reported[c("OW", "OI")] * (1 + bias)
        index <- VB * exp(sigma_index * d$index_z[k, ] - .5 * sigma_index^2)
        area_catch <- c(sum(reported[fleet_area == 1]), sum(reported[fleet_area == 2]), sum(reported[fleet_area == 3]))
        # ---- 11. 估计方法（与历史阶段相同）| estimation method (as in the history) ----
        em_enso <- if (estimation_movement == "observed") enso else 0
        Pe <- movement
        Pe[1, 2] <- Pe[1, 2] + movement_shift[1] * em_enso; Pe[1, 1] <- Pe[1, 1] - movement_shift[1] * em_enso
        Pe[2, 3] <- Pe[2, 3] + movement_shift[2] * em_enso; Pe[2, 2] <- Pe[2, 2] - movement_shift[2] * em_enso
        Nm_e <- t(Pe) %*% est
        VB_pred <- pmax(as.numeric(Nm_e %*% (sel * w)), 1e-12)
        filtered <- Nm_e * exp(index_gain * log(pmax(index, 1e-12) / VB_pred))
        SB_e <- as.numeric(filtered %*% (mat * w))
        F_e <- numeric(3); fitted <- numeric(3)
        for (a in 1:3) {
          Na <- filtered[a, ]; Ta <- area_catch[a]
          if (sum(Na * w) == 0 || Ta == 0) next
          req <- function(x) { Zx <- M + x * sel; pmin(em_Fmax, Ta / sum(sel / Zx * (-expm1(-Zx)) * Na * w)) }
          tot <- if (sum(req(em_Fmax)) - em_Fmax >= -1e-12) em_Fmax else
            uniroot(function(x) sum(req(x)) - x, c(0, em_Fmax), tol = 1e-11)$root
          F_e[a] <- req(tot)
        }
        for (a in 1:3) { Zx <- M + F_e[a] * sel; fitted[a] <- F_e[a] * sum(sel / Zx * (-expm1(-Zx)) * filtered[a, ] * w) }
        Ze <- outer(F_e, sel) + Mmat
        rec_e <- 4 * h * R0 * sum(SB_e) / (SB0 * (1 - h) + sum(SB_e) * (5 * h - 1))
        surve <- filtered * exp(-Ze)
        est <- matrix(0, 3, A); est[, 2:A] <- surve[, 1:(A - 1)]; est[, A] <- est[, A] + surve[, A]; est[, 1] <- rec_e * recruit_share
        D_est <- sum(SB_e) / SB0
        recent_catch <- rbind(recent_catch[-1, , drop = FALSE], reported)
        recent_D <- c(tail(recent_D, cycle - 1), D_est)
        # ---- 12. 本年表现指标 | this year's performance quantities ----
        D <- sum(SB) / mean(tail(U_history, unfished_window))           # 耗竭 = SB / 前 10 年未捕捞平均 | depletion
        U_history <- c(tail(U_history, unfished_window - 1), sum(SBu))
        C <- sum(catch); finite_sum <- sum(targets[controlled])
        uncapped <- !controlled & !c(FALSE, is.finite(cap), is.finite(cap), FALSE)
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
          changes = ch_WI + ch_IW, deduction = deduction, shortfall = sum(pmax(area_catch - fitted, 0)),
          catch = catch, attempted = ifelse(controlled, attempted, NA_real_))
      }
      # ---- 13. 汇总一个重复（评估期 = 第 11–20 年）| summarise one replicate (evaluation period = years 11–20) ----
      Dv <- sapply(yv, `[[`, "D"); Cv <- sapply(yv, `[[`, "C_total")
      Cf <- t(sapply(yv, `[[`, "catch")); Af <- t(sapply(yv, `[[`, "attempted"))
      e <- evaluation_years
      De <- Dv[e]; Ce <- Cv[e]
      Dt <- sapply(yv, `[[`, "D_true_fixed")[e]; Dm <- sapply(yv, `[[`, "D_est")[e]
      dec <- sapply(yv, `[[`, "decision_year")
      # 恢复：3 年滑动平均首次达到目标的年份 | recovery: first year the 3-year running mean reaches the target
      if (i0$initial_true_D >= target_D) recovered <- c(1, 0) else {
        m3 <- numeric(0)
        for (t3 in 3:n_proj) m3 <- c(m3, mean(Dv[(t3 - 2):t3]))
        hit <- which(m3 >= target_D)
        recovered <- if (length(hit)) c(1, hit[1] + 2) else c(0, n_proj)
      }
      # 配额被捕满的年份比例（只计配额为正的年份）| share of years with the quota fully caught (positive-quota years only)
      full_caught <- rep(NA_real_, 4)
      for (fi in 1:4) {
        qa <- Af[e, fi]; act <- is.finite(qa) & qa > 0
        if (any(act)) full_caught[fi] <- mean(Cf[e, fi][act] >= qa[act] * (1 - 1e-6))
      }
      rows[[length(rows) + 1]] <- data.frame(
        case_id = cell_id, mp = mp, rep = r,
        P_below_LRP = mean(De < limit_reference), P_ever_LRP = as.numeric(any(De < limit_reference)),
        P_below_42 = mean(De < diagnostic_threshold), P_ever_42 = as.numeric(any(De < diagnostic_threshold)),
        mean_D = mean(De), mean_relative_D = mean(De / i0$true_state_reference),
        mae_target = mean(abs(De - target_D)),                         # 绝对目标误差 | absolute target error
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
        relative_catch = mean(Ce) / i0$last_total_catch,              # 渔获 / 历史末年总渔获 | catch / final historical catch
        P_below_target = mean(De < target_D),
        mean_D_true_fixed = mean(Dt), mean_D_est_fixed = mean(Dm),
        estimation_bias_fixed = mean(Dm - Dt), estimation_rmse_fixed = sqrt(mean((Dm - Dt)^2)),
        mean_denominator_gap = mean(De - Dt),
        quota_fully_caught_W = full_caught[1], quota_fully_caught_OW = full_caught[2],
        quota_fully_caught_OI = full_caught[3], quota_fully_caught_E = full_caught[4],
        S1 = cell$S1, S2 = cell$S2, S3 = cell$S3, S4 = cell$S4, om_scenario = cfg$biology)
      diag_max <- pmax(diag_max, c(max(sapply(yv, `[[`, "tac_violation")), max(sapply(yv, `[[`, "capacity_violation")),
                                   max(sapply(yv, `[[`, "cap_violation")), 0, max(sapply(yv, `[[`, "shortfall"))))
      # ---- 14. 保存第 1 个重复的轨迹（图 3）| keep replicate 1's trajectory (Fig. 3) ----
      if (r == 1) traj[[length(traj) + 1]] <- data.frame(
        case_id = cell_id, mp = mp, rep = r, year = seq_len(n_proj), D = Dv,
        target_D = target_D, decision_year = dec, relative_status = sapply(yv, `[[`, "rs"),
        hcr_tac = sapply(yv, `[[`, "hcr_tac"), tac = sapply(yv, `[[`, "tac"), C_total = Cv,
        C_W = Cf[, 1], C_OW = Cf[, 2], C_OI = Cf[, 3], C_E = Cf[, 4],
        attempted_TAC_W = Af[, 1], attempted_TAC_OW = Af[, 2], attempted_TAC_OI = Af[, 3], attempted_TAC_E = Af[, 4],
        last_total_catch = i0$last_total_catch)
    }
    diag_rows[[length(diag_rows) + 1]] <- data.frame(case_id = cell_id, mp = mp,
      max_requested_TAC_violation = diag_max[["tac"]], max_capacity_violation = diag_max[["capacity"]],
      max_effort_cap_violation = diag_max[["cap"]], max_em_catch_shortfall = diag_max[["shortfall"]])
  }
  if (!VERBOSE) cat(format(Sys.time(), "%H:%M:%S"), cell_id, "done\n")
}
metrics <- bind_rows(rows); trajectories <- bind_rows(traj); diagnostics <- bind_rows(diag_rows)
# 数值检查：配额与容量约束没有被突破 | numerical checks: quota and capacity constraints hold
stopifnot(max(diagnostics$max_requested_TAC_violation) < 1e-7, max(diagnostics$max_capacity_violation) < 1e-10,
          max(diagnostics$max_effort_cap_violation) < 1e-10)

# ---- 15. 看结果 | look at the results ----
if (VERBOSE) {
  for (mp in RUN_MPS) {
    tr <- trajectories[trajectories$mp == mp, ]
    cat("\n====", mp, "| scenario", RUN_CELLS, "| replicate 1 ====\n")
    print(data.frame(year = tr$year, decision = tr$decision_year, rel_status = round(tr$relative_status, 3),
                     HCR_catch = round(tr$hcr_tac, 1), TAC = round(tr$tac, 1),
                     W = round(tr$C_W, 1), OW = round(tr$C_OW, 1), OI = round(tr$C_OI, 1), E = round(tr$C_E, 1),
                     depletion = round(tr$D, 4)), row.names = FALSE)
  }
  png(file.path(OUT, "06_depletion_rep1.png"), width = 1400, height = 900, res = 200)
  col <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7")
  plot(NA, xlim = c(1, n_proj), ylim = range(trajectories$D, trajectories$target_D), xlab = "Projection year",
       ylab = "Depletion", main = paste("Scenario", RUN_CELLS, "- replicate 1"))
  for (i in seq_along(RUN_MPS)) lines(seq_len(n_proj), trajectories$D[trajectories$mp == RUN_MPS[i]], col = col[i], lwd = 2)
  abline(h = trajectories$target_D[1], lty = 2)
  legend("bottomright", substr(RUN_MPS, 1, 4), col = col, lwd = 2, bty = "n")
  invisible(dev.off())
}
print(metrics |> group_by(case_id, mp) |>
        summarise(n = n(), depletion = mean(mean_D), relative_catch = mean(relative_catch),
                  abs_target_error = mean(mae_target), .groups = "drop"), width = Inf)
cat("06 done:", nrow(metrics), "trajectories\n")
