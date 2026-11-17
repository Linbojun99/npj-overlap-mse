# 04_initialise.R — 未捕捞平衡、可捕系数校准、参考船队份额（由参数推得的量）
# 04_initialise.R — unfished equilibrium, catchability calibration, reference fleet shares (derived quantities)
#
# 一年的种群更新顺序（全文相同）| order of one annual population step (same everywhere):
#   移动 move → 产卵潜力 spawning potential → 捕捞与自然死亡 fishing and natural mortality
#   → Beverton–Holt 补充 recruitment → 增龄 ageing (12+ 累积 accumulates)

if (!exists("sel")) source("R/03_parameters.R")
Mmat <- matrix(M, 3, A, byrow = TRUE)                  # 每区同一 M 向量 | same M vector in every region

# ---- 1. 未捕捞平衡：每年补充 R0，跑 400 年 | unfished equilibrium: recruitment R0 each year, 400 years ----
N <- matrix(0, 3, A); N[, 1] <- R0 * recruit_share
for (i in 1:400) {
  Nm <- t(movement) %*% N                              # 移动 | movement
  Z <- outer(rep(0, 3), sel) + Mmat                    # 无捕捞 | no fishing
  surv <- Nm * exp(-Z)
  N <- matrix(0, 3, A)
  N[, 2:A] <- surv[, 1:(A - 1)]                        # 增龄 | ageing
  N[, A] <- N[, A] + surv[, A]                         # 12+ 组累积 | plus group
  N[, 1] <- R0 * recruit_share                         # 补充 | recruitment
}
N0 <- N
SB0_area <- as.numeric(N0 %*% (mat * w))               # 各区未捕捞产卵潜力 | unfished spawning potential by region
SB0 <- sum(SB0_area)
cat("unfished spawning potential SB0 =", SB0, "| by region:", round(SB0_area, 1), "\n")

# ---- 2. 校准可捕系数 q：固定努力量下长期平衡耗竭 = 0.45 ----
# ---- 2. calibrate catchability q: long-run depletion at fixed effort = 0.45 ----
# 对每个候选 q，从未捕捞状态以 F = q × 容量跑 400 年，取产卵潜力 / SB0。
# For each trial q, run 400 years from unfished at F = q × capacity and take spawning potential / SB0.
steady_depletion <- function(q) {                      # uniroot 需要一个函数 | uniroot needs a function
  N <- N0
  for (i in 1:400) {
    Nm <- t(movement) %*% N
    SB <- as.numeric(Nm %*% (mat * w))
    Fy <- q * effort_end
    F_area <- c(sum(Fy[fleet_area == 1]), sum(Fy[fleet_area == 2]), sum(Fy[fleet_area == 3]))
    Z <- outer(F_area, sel) + Mmat
    rec <- 4 * h * R0 * sum(SB) / (SB0 * (1 - h) + sum(SB) * (5 * h - 1)) * 1
    surv <- Nm * exp(-Z)
    N <- matrix(0, 3, A); N[, 2:A] <- surv[, 1:(A - 1)]; N[, A] <- N[, A] + surv[, A]
    N[, 1] <- rec * recruit_share
  }
  N
}
q_end <- uniroot(function(q) sum(as.numeric(steady_depletion(q) %*% (mat * w))) / SB0 - start_depletion,
                 c(.001, 3), tol = 1e-10)$root
q_initial <- q_end / effort_creep^(n_hist - 1)          # 历史第 1 年的 q | q in the first history year
# 历史期努力量每年除以 1.01，抵消 q 增长，使历史 F 不变 | history effort offsets catchability growth, F constant
history_effort <- outer(effort_creep^(n_hist - seq_len(n_hist)), effort_end)
colnames(history_effort) <- fleets
N_start <- steady_depletion(q_end)                     # 历史第 1 年的种群 | population at the start of the history
cat("q_end =", q_end, "| q_initial =", q_initial, "\n")

# ---- 3. 参考船队份额：平衡耗竭 0.56 时各船队渔获占比（MP-A、C、D 的分配比例）----
# ---- 3. reference fleet shares: catch shares at equilibrium depletion 0.56 (allocation of MP-A, C, D) ----
equilibrium_catch <- function(scale) {                 # 300 年后的产卵潜力与各船队渔获 | spawning potential and catch after 300 years
  N <- N0
  for (i in 1:300) {
    Nm <- t(movement) %*% N
    SB <- as.numeric(Nm %*% (mat * w))
    Fy <- q_end * effort_end * scale
    F_area <- c(sum(Fy[fleet_area == 1]), sum(Fy[fleet_area == 2]), sum(Fy[fleet_area == 3]))
    Z <- outer(F_area, sel) + Mmat
    g <- c(sum(sel / Z[1, ] * (-expm1(-Z[1, ])) * Nm[1, ] * w),     # 每单位 F 的渔获 | catch per unit F
           sum(sel / Z[2, ] * (-expm1(-Z[2, ])) * Nm[2, ] * w),
           sum(sel / Z[3, ] * (-expm1(-Z[3, ])) * Nm[3, ] * w))
    catch <- Fy * g[fleet_area]                        # Baranov 渔获 | Baranov catch
    rec <- 4 * h * R0 * sum(SB) / (SB0 * (1 - h) + sum(SB) * (5 * h - 1)) * 1
    surv <- Nm * exp(-Z)
    N <- matrix(0, 3, A); N[, 2:A] <- surv[, 1:(A - 1)]; N[, A] <- N[, A] + surv[, A]
    N[, 1] <- rec * recruit_share
  }
  list(SB = SB, catch = setNames(catch, fleets))
}
scale <- uniroot(function(x) sum(equilibrium_catch(x)$SB) / SB0 - illustrative_target, c(.001, 3), tol = 1e-8)$root
eq <- equilibrium_catch(scale)
C_target <- sum(eq$catch)
share_ref <- eq$catch / C_target                       # W, OW, OI, E 的参考份额 | reference shares
managed_share <- sum(share_ref[c("W", "OW")])          # MP 管理的份额（W+OW）| share managed through the MP
overlap_effort_cap <- overlap_cap_fraction * sum(effort_end[c("OW", "OI")])   # MP-D 上限 | MP-D cap
print(round(share_ref, 4))
cat("MP-D overlap effort cap =", overlap_effort_cap, "\n")

cat("04 done\n")
