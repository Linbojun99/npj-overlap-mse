# 11_external_catch.R — 补充实验：MP 范围外渔获，两侧对称（补充说明 7、补充图 5、补充表 14）
# 11_external_catch.R — supplementary experiment: catch outside the MP's scope on both sides
#                       (Supplementary Note 7, Supplementary Fig. 5, Supplementary Table 14)
#
# WCPFC 的官方评估（WCPFC22-2025-21）对两支不受 MP 管理的渔业设了固定渔获：东太（重疊区以外，船队 E）
# 与中西太公约区赤道至 10°S 的延绳钓（这里加的船队 T，R1，历史渔获 = R1+R2 的 12%）。
# 两支船队的容量各取每年 −3%、不变、+3%，3 × 3 全交叉 × S1–S3 的 12 个组合 × 4 个 MP。
# The WCPFC evaluation (WCPFC22-2025-21) holds two fisheries outside the MP at fixed catches: the EPO outside
# the overlap area (fleet E) and the WCPFC longline fishery between the equator and 10°S (fleet T added here,
# in R1, 12% of the historical R1+R2 catch). Each fleet's capacity changes by −3%, 0 or +3% per year;
# 3 × 3 levels × 12 combinations of S1–S3 × 4 MPs.
#   N_REP=20 Rscript R/11_external_catch.R     # 约 6 分钟 | about 6 min
#   N_REP=300 Rscript R/11_external_catch.R    # 论文设置（约 1.5 小时）| as in the Article (about 1.5 h)

if (!exists("N_REP")) N_REP <- as.integer(Sys.getenv("N_REP", "20"))
CONFIG <- "reference"
source("R/00_settings.R")
source("R/03_parameters.R"); source("R/04_initialise.R")

# ---- 1. 新增船队 T 的参数 | parameters of the added fleet T ----
t_share <- .12        # WCPFC22-2025-21 第 3 页：赤道至 10°S 约占中西太近年渔获 12% | p. 3: about 12% of recent WCPFC-CA catch
growth <- .03         # 与 S4 同幅度，两支船队、两个方向相同 | same amplitude as S4, both fleets, both directions
xfleets <- c("W", "OW", "OI", "E", "T")
xarea <- c(1L, 2L, 2L, 3L, 1L)                    # T 在 R1 与 W 共用种群 | T fishes R1 with W
# T 从 W 容量中划出 ext_v，使历史 T 渔获 = 12% × (W+OW+OI)；R1 的历史 F 不变，历史种群与主实验相同。
# T takes a fraction ext_v of W's capacity so that historical T catch = 12% × (W+OW+OI); R1's historical F is
# unchanged, so the historical population equals the main experiment's.
Nm <- t(movement) %*% N_start                     # 历史末年的平衡状态 | fished steady state
Fy <- q_end * effort_end
F_area <- c(sum(Fy[fleet_area == 1]), sum(Fy[fleet_area == 2]), sum(Fy[fleet_area == 3]))
Z <- outer(F_area, sel) + Mmat
g <- c(sum(sel / Z[1, ] * (-expm1(-Z[1, ])) * Nm[1, ] * w), sum(sel / Z[2, ] * (-expm1(-Z[2, ])) * Nm[2, ] * w),
       sum(sel / Z[3, ] * (-expm1(-Z[3, ])) * Nm[3, ] * w))
c0 <- setNames(Fy * g[fleet_area], fleets)
ext_v <- unname(t_share * sum(c0[c("W", "OW", "OI")]) / c0[["W"]])
# 四船队向量 → 五船队：W 拆成 W(1−v) 与 T = W·v | four-fleet vector → five fleets
xshare_ref <- setNames(c(share_ref[[1]] * (1 - ext_v), as.numeric(share_ref[2:4]), share_ref[[1]] * ext_v), xfleets)
xmanaged <- sum(xshare_ref[c("W", "OW")])         # MP 管理的份额（不含 T）| share managed through the MP (T excluded)
cat("fraction of W capacity assigned to T:", round(ext_v, 4), "\n"); print(round(xshare_ref, 4))

# ---- 2. 随机数与历史（同 05，历史渔获拆成五船队）| random numbers and history (as step 05, catch split into five fleets) ----
source("R/05_history.R")
for (r in seq_len(N_REP)) {
  h4 <- init[[r]]$catch_all                                        # 30 年 × 4 船队观测渔获 | 30 years × 4 fleets
  ch <- cbind(W = h4[, "W"] * (1 - ext_v), h4[, c("OW", "OI", "E")], T = h4[, "W"] * ext_v)   # W 拆出 T | split T out of W
  ci <- match(catch_reference_years, history_years)
  init[[r]]$catch_reference <- mean(rowSums(ch[ci, c("W", "OW"), drop = FALSE]))   # 只计 W、OW | W and OW only
  init[[r]]$change_base_catch <- sum(ch[n_hist, c("W", "OW")])
  init[[r]]$last_total_catch <- sum(ch[n_hist, ])
  init[[r]]$catch_history <- tail(ch, cycle)
}

# ---- 3. 闭环（与 06 相同，船队变为五支）| closed loop (as step 06, with five fleets) ----
x1 <- hcr_points[1, 1]; y1 <- hcr_points[1, 2]; x2 <- hcr_points[2, 1]; y2 <- hcr_points[2, 2]
x3 <- hcr_points[3, 1]; x4 <- hcr_points[4, 1]; y4 <- hcr_points[4, 2]
cells <- design[design$S4 == 0, ]                 # 12 个组合（S4 由 E 的水平取代）| 12 cells (S4 replaced by E levels)
rows <- list()
for (lev_T in -1:1) for (lev_E in -1:1) {         # E 变化最快 | E varies fastest
  for (gi in seq_len(nrow(cells))) {
    cell <- cells[gi, ]; S1 <- cell$S1 == 1; S2 <- cell$S2; S3 <- cell$S3 == 1
    bias_centre <- abs(misreporting) * S2
    for (mp in mps) {
      mp_i <- match(mp, mps)
      for (r in seq_len(N_REP)) {
        d <- draws[[r]]; i0 <- init[[r]]
        N <- i0$N; U <- i0$U; est <- i0$est; U_history <- i0$U_history
        target_D <- target_relative * i0$true_state_reference
        recent_catch <- i0$catch_history; recent_D <- i0$depletion_history
        nW <- round(n_segments * effort_end[["OW"]] / sum(effort_end[c("OW", "OI")]))
        seg_measure <- c(rep("W", nW), rep("I", n_segments - nW))
        seg_capacity <- rep(sum(effort_end[c("OW", "OI")]) / n_segments, n_segments)
        seg_last <- rep(-min_measure_years, n_segments)
        external_capacity <- effort_end[["E"]]
        tropical_capacity <- ext_v * effort_end[["W"]]
        tac <- i0$change_base_catch; first_tac <- tac
        targets <- setNames(rep(Inf, 5), xfleets)
        enso <- 0; deduction <- 0
        yv <- list()
        for (y in seq_len(n_proj)) {
          k <- n_hist + y; q <- q_end * effort_creep^y
          if (S3) enso <- c(-1, 0, 1)[which(d$enso_u[y] <= cumsum(env_transition[match(enso, -1:1), ]))[1]]
          # E 与 T 的容量按各自水平变化 | E and T capacities change by their levels
          external_capacity <- external_capacity * (1 + lev_E * growth)
          tropical_capacity <- tropical_capacity * (1 + lev_T * growth)
          capacity <- c(W = effort_end[["W"]] * (1 - ext_v), OW = sum(seg_capacity[seg_measure == "W"]),
                        OI = sum(seg_capacity[seg_measure == "I"]), E = external_capacity, T = tropical_capacity)
          cap <- if (mp_i == 4) overlap_effort_cap else Inf
          ch_n <- 0L
          decision_year <- (y - 1) %% cycle == 0
          if (decision_year) {
            decision_D <- mean(recent_D); rs <- decision_D / i0$state_reference
            zz <- pmin(1, pmax(0, (rs - x1) / (x2 - x1)))
            rise <- -expm1(-hcr_curvature * (x2 - x1) * zz) / -expm1(-hcr_curvature * (x2 - x1))
            hcr_tac <- i0$catch_reference * ifelse(rs < x2, y1 + (y2 - y1) * rise,
                                                   y2 + (y4 - y2) * pmin(1, pmax(0, (rs - x3) / (x4 - x3))))
            previous_tac <- tac
            lower_lim <- change_bounds[1] * previous_tac; upper_lim <- change_bounds[2] * previous_tac
            tac <- min(upper_lim, max(lower_lim, hcr_tac)); cl_lower <- hcr_tac < lower_lim; cl_upper <- hcr_tac > upper_lim
            # 分配：MP-B 的共同 TAC 也覆盖 T；MP-C 的扣减计入 OI、E、T | MP-B's joint TAC covers T; MP-C deducts OI, E and T
            targets <- setNames(rep(Inf, 5), xfleets); deduction <- 0
            if (mp_i == 2) {
              obs_sh <- colSums(recent_catch)
              obs_sh <- if (sum(obs_sh) > 0) obs_sh / sum(obs_sh) else xshare_ref
              targets <- tac / xmanaged * obs_sh
            } else {
              targets[c("W", "OW")] <- tac * xshare_ref[c("W", "OW")] / xmanaged
              if (mp_i == 3) {
                external <- mean(rowSums(recent_catch[, c("OI", "E", "T"), drop = FALSE]))
                deduction <- max(0, external - previous_tac / xmanaged * sum(xshare_ref[c("OI", "E", "T")]))
                total <- sum(targets[c("W", "OW")]); deduction <- min(deduction, total)
                targets[c("W", "OW")] <- if (total > 0) targets[c("W", "OW")] * (1 - deduction / total) else 0
              }
            }
            if (mp_i == 2) { tot <- sum(targets[c("OW", "OI")]); den <- sum(capacity[c("OW", "OI")])
                             targets[c("OW", "OI")] <- if (den > 0) tot * capacity[c("OW", "OI")] / den else 0 }
            Nm_p <- t(movement) %*% est
            available <- capacity; ov <- sum(available[c("OW", "OI")])
            if (ov > cap) available[c("OW", "OI")] <- available[c("OW", "OI")] * cap / ov
            Fp <- setNames(numeric(5), xfleets)
            for (a in 1:3) {   # R1 中 W 与 T 联合求解 | W and T solved jointly in R1
              j <- which(xarea == a); Fcap <- q * available[j]; Tj <- targets[j]; Na <- Nm_p[a, ]
              if (sum(Fcap) == 0 || sum(Na * w) == 0 || all(Tj == 0)) { Fp[j] <- Fcap * 0; next }
              req <- function(x) { Zx <- M + x * sel; pmin(Fcap, Tj / sum(sel / Zx * (-expm1(-Zx)) * Na * w)) }
              tot <- if (sum(req(sum(Fcap))) - sum(Fcap) >= -1e-12) sum(Fcap) else
                uniroot(function(x) sum(req(x)) - x, c(0, sum(Fcap)), tol = 1e-11)$root
              Fp[j] <- req(tot)
            }
            pred_effort <- Fp / q
            use_OW <- if (capacity[["OW"]] > 0) pred_effort[["OW"]] / capacity[["OW"]] else as.numeric(targets[["OW"]] > 0)
            use_OI <- if (capacity[["OI"]] > 0) pred_effort[["OI"]] / capacity[["OI"]] else as.numeric(targets[["OI"]] > 0)
            gap <- use_OI - use_OW
            p_WI <- if (gap <= switch_cost) 0 else plogis(switch_intercept + switch_slope * (gap - switch_cost))
            p_IW <- if (-gap <= switch_cost) 0 else plogis(switch_intercept + switch_slope * (-gap - switch_cost))
            if (S1) {
              eligible <- y - seg_last >= min_measure_years
              ch <- eligible & d$measure_change_u[y, ] < switch_strength * ifelse(seg_measure == "W", p_WI, p_IW)
              old <- seg_measure
              if (any(ch)) { seg_measure[ch] <- ifelse(old[ch] == "W", "I", "W"); seg_last[ch] <- as.integer(y) }
              ch_n <- sum(ch & old == "W") + sum(ch & old == "I")
              capacity <- c(W = effort_end[["W"]] * (1 - ext_v), OW = sum(seg_capacity[seg_measure == "W"]),
                            OI = sum(seg_capacity[seg_measure == "I"]), E = external_capacity, T = tropical_capacity)
              if (mp_i == 2) { tot <- sum(targets[c("OW", "OI")]); den <- sum(capacity[c("OW", "OI")])
                               targets[c("OW", "OI")] <- if (den > 0) tot * capacity[c("OW", "OI")] / den else 0 }
            }
          }
          # 实施：T 与 W 共用实施误差（参考配置为 0）| implementation: T shares W's error (zero here)
          impl <- exp(sigma_impl * d$implementation_z[y, c(1:4, 1L)] - .5 * sigma_impl^2)
          controlled <- is.finite(targets)
          attempted <- targets; attempted[controlled] <- targets[controlled] * impl[controlled]
          Py <- movement
          Py[1, 2] <- Py[1, 2] + movement_shift[1] * enso; Py[1, 1] <- Py[1, 1] - movement_shift[1] * enso
          Py[2, 3] <- Py[2, 3] + movement_shift[2] * enso; Py[2, 2] <- Py[2, 2] - movement_shift[2] * enso
          Nm <- t(Py) %*% N
          available <- capacity; ov <- sum(available[c("OW", "OI")])
          if (ov > cap) available[c("OW", "OI")] <- available[c("OW", "OI")] * cap / ov
          Fy <- setNames(numeric(5), xfleets)
          for (a in 1:3) {
            j <- which(xarea == a); Fcap <- q * available[j]; Tj <- attempted[j]; Na <- Nm[a, ]
            if (sum(Fcap) == 0 || sum(Na * w) == 0 || all(Tj == 0)) { Fy[j] <- Fcap * 0; next }
            req <- function(x) { Zx <- M + x * sel; pmin(Fcap, Tj / sum(sel / Zx * (-expm1(-Zx)) * Na * w)) }
            tot <- if (sum(req(sum(Fcap))) - sum(Fcap) >= -1e-12) sum(Fcap) else
              uniroot(function(x) sum(req(x)) - x, c(0, sum(Fcap)), tol = 1e-11)$root
            Fy[j] <- req(tot)
          }
          effort <- Fy / q
          # 真实种群 | true population
          rmult <- exp(sigma_R * d$recruitment_z[k] - .5 * sigma_R^2)
          SB <- as.numeric(Nm %*% (mat * w)); VB <- as.numeric(Nm %*% (sel * w))
          F_area <- c(sum(Fy[xarea == 1]), sum(Fy[xarea == 2]), sum(Fy[xarea == 3]))
          Z <- outer(F_area, sel) + Mmat
          g <- c(sum(sel / Z[1, ] * (-expm1(-Z[1, ])) * Nm[1, ] * w), sum(sel / Z[2, ] * (-expm1(-Z[2, ])) * Nm[2, ] * w),
                 sum(sel / Z[3, ] * (-expm1(-Z[3, ])) * Nm[3, ] * w))
          catch <- setNames(Fy * g[xarea], xfleets)
          rec <- 4 * h * R0 * sum(SB) / (SB0 * (1 - h) + sum(SB) * (5 * h - 1)) * rmult
          surv <- Nm * exp(-Z)
          N <- matrix(0, 3, A); N[, 2:A] <- surv[, 1:(A - 1)]; N[, A] <- N[, A] + surv[, A]; N[, 1] <- rec * recruit_share
          Um <- t(Py) %*% U; SBu <- as.numeric(Um %*% (mat * w))
          recu <- 4 * h * R0 * sum(SBu) / (SB0 * (1 - h) + sum(SBu) * (5 * h - 1)) * rmult
          survu <- Um * exp(-(outer(c(0, 0, 0), sel) + Mmat))
          U <- matrix(0, 3, A); U[, 2:A] <- survu[, 1:(A - 1)]; U[, A] <- U[, A] + survu[, A]; U[, 1] <- recu * recruit_share
          # 观测：T 与 W 一同上报 WCPFC，再按真实渔获比例拆开 | T reports with W to WCPFC, then split by true catch
          parent <- setNames(c(catch[["W"]] + catch[["T"]], as.numeric(catch[2:4])), fleets)
          bias <- if (S2 != 0) bias_centre * (.5 + d$reporting_u[k]) else 0
          reported <- parent * exp(sigma_catch * d$catch_z[k, ] - .5 * sigma_catch^2)
          reported[c("OW", "OI")] <- reported[c("OW", "OI")] * (1 + bias)
          index <- VB * exp(sigma_index * d$index_z[k, ] - .5 * sigma_index^2)
          area_catch <- c(sum(reported[fleet_area == 1]), sum(reported[fleet_area == 2]), sum(reported[fleet_area == 3]))
          tf <- if (parent[["W"]] > 0) catch[["T"]] / parent[["W"]] else 0
          reported5 <- setNames(c(reported[[1]] * (1 - tf), as.numeric(reported[2:4]), reported[[1]] * tf), xfleets)
          # 估计方法（与 06 相同，只看分区渔获）| estimation method (as step 06; regional catches only)
          Nm_e <- t(movement) %*% est
          VB_pred <- pmax(as.numeric(Nm_e %*% (sel * w)), 1e-12)
          filtered <- Nm_e * exp(index_gain * log(pmax(index, 1e-12) / VB_pred))
          SB_e <- as.numeric(filtered %*% (mat * w))
          F_e <- numeric(3)
          for (a in 1:3) {
            Na <- filtered[a, ]; Ta <- area_catch[a]
            if (sum(Na * w) == 0 || Ta == 0) next
            req <- function(x) { Zx <- M + x * sel; pmin(em_Fmax, Ta / sum(sel / Zx * (-expm1(-Zx)) * Na * w)) }
            tot <- if (sum(req(em_Fmax)) - em_Fmax >= -1e-12) em_Fmax else uniroot(function(x) sum(req(x)) - x, c(0, em_Fmax), tol = 1e-11)$root
            F_e[a] <- req(tot)
          }
          surve <- filtered * exp(-(outer(F_e, sel) + Mmat))
          rec_e <- 4 * h * R0 * sum(SB_e) / (SB0 * (1 - h) + sum(SB_e) * (5 * h - 1))
          est <- matrix(0, 3, A); est[, 2:A] <- surve[, 1:(A - 1)]; est[, A] <- est[, A] + surve[, A]; est[, 1] <- rec_e * recruit_share
          D_est <- sum(SB_e) / SB0
          recent_catch <- rbind(recent_catch[-1, , drop = FALSE], reported5)
          recent_D <- c(tail(recent_D, cycle - 1), D_est)
          D <- sum(SB) / mean(tail(U_history, unfished_window))
          U_history <- c(tail(U_history, unfished_window - 1), sum(SBu))
          C <- sum(catch); finite_sum <- sum(targets[controlled])
          uncapped <- !controlled & !c(FALSE, is.finite(cap), is.finite(cap), FALSE, FALSE)
          yv[[y]] <- list(D = D, D_true_fixed = sum(SB) / SB0, D_est = D_est, D_overlap = SB[2] / SBu[2],
            decision_year = decision_year, cl_lower = cl_lower, cl_upper = cl_upper, tac = tac, C_total = C,
            not_managed = if (C > 0) sum(catch[uncapped]) / C else 0,
            outside_W = if (C > 0) sum(catch[c("OI", "E")]) / C else 0,
            exceed = if (finite_sum > 0) sum(pmax(catch[controlled] - targets[controlled], 0)) / finite_sum else 0,
            uncaught = if (finite_sum > 0) sum(pmax(targets[controlled] - catch[controlled], 0)) / finite_sum else 0,
            E_overlap = sum(effort[c("OW", "OI")]),
            cap_violation = if (is.finite(cap)) max(0, sum(effort[c("OW", "OI")]) - cap) else 0,
            changes = ch_n, deduction = deduction, catch = catch, attempted = ifelse(controlled, attempted, NA_real_))
        }
        # 汇总（同 06，加上 T 与按管理范围的渔获）| summary (as step 06, plus T and catch by scope)
        Dv <- sapply(yv, `[[`, "D"); Cv <- sapply(yv, `[[`, "C_total")
        Cf <- t(sapply(yv, `[[`, "catch")); Af <- t(sapply(yv, `[[`, "attempted"))
        e <- evaluation_years; De <- Dv[e]; Ce <- Cv[e]
        Dt <- sapply(yv, `[[`, "D_true_fixed")[e]; Dm <- sapply(yv, `[[`, "D_est")[e]; dec <- sapply(yv, `[[`, "decision_year")
        if (i0$initial_true_D >= target_D) recovered <- c(1, 0) else {   # 3 年滑动平均达到目标 | 3-year running mean reaches target
          m3 <- numeric(0)
          for (t3 in 3:n_proj) m3 <- c(m3, mean(Dv[(t3 - 2):t3]))
          hit <- which(m3 >= target_D)
          recovered <- if (length(hit)) c(1, hit[1] + 2) else c(0, n_proj)
        }
        full_caught <- rep(NA_real_, 4)                              # 配额被捕满的年份比例 | share of years quota fully caught
        for (fi in 1:4) {
          qa <- Af[e, fi]; act <- is.finite(qa) & qa > 0
          if (any(act)) full_caught[fi] <- mean(Cf[e, fi][act] >= qa[act] * (1 - 1e-6))
        }
        anchor <- i0$last_total_catch
        rows[[length(rows) + 1]] <- data.frame(
          case_id = cell$case_id, mp = mp, rep = r,
          P_below_LRP = mean(De < limit_reference), P_ever_LRP = as.numeric(any(De < limit_reference)),
          P_below_42 = mean(De < diagnostic_threshold), P_ever_42 = as.numeric(any(De < diagnostic_threshold)),
          mean_D = mean(De), mean_relative_D = mean(De / i0$true_state_reference), mae_target = mean(abs(De - target_D)),
          mean_target = mean(rep(target_D, length(e))), signed_target_gap = mean(De - target_D), mean_catch = mean(Ce),
          AAV = if (sum(head(Ce, -1)) > 0) sum(abs(diff(Ce))) / sum(head(Ce, -1)) else NA_real_,
          not_managed_fraction = mean(sapply(yv, `[[`, "not_managed")[e]), outside_W_fraction = mean(sapply(yv, `[[`, "outside_W")[e]),
          D_overlap = mean(sapply(yv, `[[`, "D_overlap")[e]), overlap_effort = mean(sapply(yv, `[[`, "E_overlap")[e]),
          quota_exceedance_ratio = mean(sapply(yv, `[[`, "exceed")[e]), uncaught_quota_ratio = mean(sapply(yv, `[[`, "uncaught")[e]),
          effort_cap_exceedance = mean(sapply(yv, `[[`, "cap_violation")[e] > 1e-8),
          measure_changes_in_window = sum(sapply(yv, `[[`, "changes")[e]), mean_deduction = mean(sapply(yv, `[[`, "deduction")[e]),
          catch_W = mean(Cf[e, 1]), catch_OW = mean(Cf[e, 2]), catch_OI = mean(Cf[e, 3]), catch_E = mean(Cf[e, 4]),
          change_limit_lower_frequency = mean(sapply(yv, `[[`, "cl_lower")[dec]),
          change_limit_upper_frequency = mean(sapply(yv, `[[`, "cl_upper")[dec]),
          change_base_initial_cut = max(0, 1 - first_tac / i0$change_base_catch),
          first_tac_cut = max(0, 1 - yv[[1]]$tac / i0$change_base_catch),
          early_catch_reduction = 1 - mean(Cv[1:3]) / anchor,
          recovery_event = recovered[1], recovery_time_restricted = recovered[2],
          relative_catch = mean(Ce) / anchor, P_below_target = mean(De < target_D),
          mean_D_true_fixed = mean(Dt), mean_D_est_fixed = mean(Dm),
          estimation_bias_fixed = mean(Dm - Dt), estimation_rmse_fixed = sqrt(mean((Dm - Dt)^2)), mean_denominator_gap = mean(De - Dt),
          quota_fully_caught_W = full_caught[1], quota_fully_caught_OW = full_caught[2],
          quota_fully_caught_OI = full_caught[3], quota_fully_caught_E = full_caught[4],
          catch_T = mean(Cf[e, 5]), catch_MP_fleets = mean(Cf[e, 1] + Cf[e, 2]),               # MP 管理的 W+OW | W+OW managed by the MP
          catch_outside_MP = mean(Cf[e, 3] + Cf[e, 4] + Cf[e, 5]),                              # 范围外 OI+E+T | outside the MP
          relative_catch_T = mean(Cf[e, 5]) / anchor, relative_catch_MP_fleets = mean(Cf[e, 1] + Cf[e, 2]) / anchor,
          relative_catch_outside_MP = mean(Cf[e, 3] + Cf[e, 4] + Cf[e, 5]) / anchor,
          S1 = cell$S1, S2 = cell$S2, S3 = cell$S3, E_level = lev_E, T_level = lev_T)
      }
    }
  }
  cat(format(Sys.time(), "%H:%M:%S"), "E level", lev_E, "| T level", lev_T, "done\n")
}
metrics <- bind_rows(rows)

# ---- 4. 保存 | save ----
dest <- file.path(OUT, "11_external_catch"); dir.create(dest, showWarnings = FALSE)
con <- gzfile(file.path(dest, "replicate_metrics.csv.gz"), "wt"); write.csv(metrics, con, row.names = FALSE, na = ""); close(con)

# ---- 5. 效应：相对两支船队都不变，同一组合、MP、重复内配对（补充表 14、补充图 5）----
# ---- 5. effects relative to both fleets constant, paired within cell, MP and replicate (Supplementary Table 14, Fig. 5) ----
ref <- metrics |> filter(E_level == 0, T_level == 0) |>
  select(case_id, mp, rep, D0 = mean_D, MP0 = relative_catch_MP_fleets, C0 = relative_catch, AE0 = mae_target)
eff <- metrics |> filter(E_level != 0 | T_level != 0) |> inner_join(ref, by = c("case_id", "mp", "rep")) |>
  transmute(mp = substr(mp, 1, 4), E_level, T_level, rep,
            mean_D = mean_D - D0, relative_catch_MP_fleets = relative_catch_MP_fleets - MP0,
            relative_catch = relative_catch - C0, mae_target = mae_target - AE0) |>
  pivot_longer(c(mean_D, relative_catch_MP_fleets, relative_catch, mae_target), names_to = "metric") |>
  group_by(mp, E_level, T_level, metric, rep) |> summarise(value = mean(value), .groups = "drop") |>   # 12 个组合内平均 | over the 12 cells
  group_by(mp, E_level, T_level, metric) |>
  summarise(effect = mean(value), mcse = sd(value) / sqrt(n()), n_rep = n(), .groups = "drop")
write.csv(eff, file.path(dest, "marginal_effects.csv"), row.names = FALSE)
cat("\nchange in depletion (percentage points):\n")
print(as.data.frame(eff |> filter(metric == "mean_D", E_level == 0 | T_level == 0 | E_level == T_level) |>
  mutate(points = round(100 * effect, 2), level = paste0("E", E_level, " T", T_level)) |>
  select(level, mp, points) |> pivot_wider(names_from = mp, values_from = points)))
cat("trajectories ever below the 20% limit reference point:", sum(metrics$P_ever_LRP), "of", nrow(metrics), "\n")
cat("11 done:", dest, "\n")
