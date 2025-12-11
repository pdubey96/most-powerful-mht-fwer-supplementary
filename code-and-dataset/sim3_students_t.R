# ==============================================================================
# Script Name: sim3_students_t.R
# Purpose:     Simulation Study 3: Heavy-Tailed Alternatives (Student's t).
#              Tests the robustness of Algorithm 1 when the deviation from the
#              null is driven by tail weight (df) rather than a mean shift.
#              Reproduces the results for Section 5.3 of the manuscript.
#
# Simulation Setting:
#   - Hypotheses: K=3 independent tests.
#   - H0: X_k ~ N(0, 1).
#   - HA: X_k ~ t_df (Student's t with degrees of freedom df).
#   - p-values: Two-sided test p = 2 * pnorm(-abs(X)).
#
# Methodology:
#   1. Vary df from 2 (heavy tails) to 20 (near normal).
#   2. Algorithm 1: Compute optimal dual variables (mu) for each df to
#      maximize Average Power (Pi_3) subject to FWER <= 0.05.
#      - Uses specific g(u) ratio: dt(x, df) / dnorm(x).
#   3. Monte Carlo: Simulate 120,000 replications to estimate power.
#   4. Comparison: Benchmark against Bonferroni, Holm, Hochberg, Hommel,
#      and Romano-Wolf (stepdown minP).
#
# Output:      - Table A: Average Power (Pi_3).
#              - Table B: Minimal Power (Pi_any).
#              - PDF Plot: 'sim_v8.1.pdf' visualizing the power curves.
# ==============================================================================

set.seed(2025)

## -----------------------------
## Global settings
## -----------------------------
K         <- 3
alpha     <- 0.05
dfs       <- seq(2, 20, by = 2) # Vary degrees of freedom
nrep      <- 120000    # Monte Carlo reps for power
Ngrid     <- 120000    # MC grid size for integrals on Q
tol_mu    <- 5e-5      # Tolerance for coordinate update
outer_tol <- 5e-5      # Tolerance for outer loop convergence
max_outer <- 60        # Max iterations

## Bracketing / stability controls for root-finding
MU_LOWER0 <- 0.0
SCAN_K    <- 40
SCAN_BASE <- 2.0
MU_CAP    <- 1e6
DBG       <- FALSE

## -----------------------------
## t-distribution helpers
## -----------------------------
Phi     <- pnorm
PhiInv  <- qnorm

## Alt p-value density g(u)
## Derivation:
##   Let x_p = qnorm(1 - p/2) be the critical value under H0 (Normal).
##   The likelihood ratio is f_A(x_p) / f_0(x_p).
##   f_A is the t-density with df degrees of freedom.
##   f_0 is the standard normal density.
g_from_u_t <- function(u, df) {
  # Clamp u to avoid Inf at u=0
  u_clamped <- pmax(pmin(u, 1 - 1e-16), 1e-16)
  xp <- PhiInv(1 - u_clamped / 2)
  
  f_A <- dt(xp, df)
  f_0 <- dnorm(xp, 0, 1)
  
  g_val <- f_A / f_0
  g_val
}

## t-distribution sampler
## Used for Monte Carlo power estimation under HA.
rg_t <- function(n, df) {
  rt(n, df)
}


## -----------------------------
## Utilities
## -----------------------------
safe_gt0 <- function(v) { w <- v > 0; w[is.na(w)] <- FALSE; w }

## Monte Carlo grid on Q (ordered simplex: u1 < u2 < u3)
sample_Q <- function(N, K=3) {
  m <- matrix(runif(N * K), ncol = K)
  t(apply(m, 1, sort))
}

## -----------------------------
## Per-df cache on shared Q grid
## -----------------------------
## Pre-calculates g(u) and the coefficients a_i, b_{l,i} for the optimization
make_df_cache <- function(df, Qgrid) {
  u1 <- Qgrid[,1]; u2 <- Qgrid[,2]; u3 <- Qgrid[,3]
  g1 <- g_from_u_t(u1, df)
  g2 <- g_from_u_t(u2, df)
  g3 <- g_from_u_t(u3, df)
  
  ## Objective Coefficients a_i(u) for Average Power (Pi_3)
  a  <- 2 * (g1 * g2 * g3)
  
  ## Constraint Coefficients b_{l,i}(u) for FWER control (from Lemma 3.1)
  g2p3 <- 2 * (g2 + g3)    # Coeff for FWER1 constraint (D1)
  g1_  <- 2 * g1          # Coeff for FWER1 constraint (D2)
  g23  <- 2 * (g2 * g3)    # Coeff for FWER2 constraint (D1)
  g13  <- 2 * (g1 * g3)    # Coeff for FWER2 constraint (D2)
  g12  <- 2 * (g1 * g2)    # Coeff for FWER2 constraint (D3)
  
  list(u1=u1,u2=u2,u3=u3,
       g1=g1,g2=g2,g3=g3,
       a=a, g2p3=g2p3, g1_=g1_, g23=g23, g13=g13, g12=g12)
}

## -----------------------------
## R_i for Π3 objective
## -----------------------------
## Computes the Lagrangian components R_i(mu, u) = a_i - sum(mu * b_i)
Ri_all_Pi3 <- function(mu, cache) {
  mu0 <- mu[1]; mu1 <- mu[2]; mu2 <- mu[3]
  R1 <- cache$a - 6*mu0                 - cache$g2p3*mu1 - cache$g23*mu2
  R2 <- cache$a                         - cache$g1_ *mu1 - cache$g13*mu2
  R3 <- cache$a                                          - cache$g12*mu2
  list(R1=R1, R2=R2, R3=R3)
}

## Determines the optimal decision policy D^mu based on R_i values
alphas_and_D <- function(R) {
  S1 <- R$R1
  S2 <- R$R1 + R$R2
  S3 <- R$R1 + R$R2 + R$R3
  a1 <- safe_gt0(S1) | safe_gt0(S2) | safe_gt0(S3)  # α1 indicator
  a2 <- safe_gt0(R$R2) | safe_gt0(R$R2 + R$R3)      # α2 indicator
  a3 <- safe_gt0(R$R3)                             # α3 indicator
  D1 <- a1
  D2 <- a1 & a2
  D3 <- a1 & a2 & a3
  list(D1=D1, D2=D2, D3=D3)
}

## -----------------------------
## TRUE FWER values at μ
## -----------------------------
## Evaluates the constraint integrals to check for slackness/equality
fwer_values_from_R <- function(R, cache) {
  D <- alphas_and_D(R)
  g1 <- cache$g1; g2 <- cache$g2; g3 <- cache$g3
  g1[!is.finite(g1)] <- 0; g2[!is.finite(g2)] <- 0; g3[!is.finite(g3)] <- 0
  
  # Integrals over Q are approximated by means over the sample grid
  F0 <- mean(D$D1)
  F1 <- (1/3) * mean( D$D1 * (g2 + g3) + D$D2 * g1 )
  F2 <- (1/3) * mean( D$D1 * (g2 * g3) + D$D2 * (g1 * g3) + D$D3 * (g1 * g2) )
  c(FWER0=F0, FWER1=F1, FWER2=F2)
}


## -----------------------------
## Algorithm 2: ComputeCoordinateMu
## -----------------------------
## Implements the Bracketing + Bisection strategy for 1D root finding
ComputeCoordinateMu <- function(Fgamma, muA_fixed, muB_fixed, alpha,
                                delta = tol_mu, MaxIter_b = 80,
                                U_s = 1e-8, U_f = SCAN_BASE, U_max = MU_CAP) {
  L <- 0.0
  mu_coord <- 0.0
  flag <- 0L
  msg  <- ""
  
  F0 <- Fgamma(L, muA_fixed, muB_fixed)
  
  if (is.finite(F0) && abs(F0 - alpha) < delta) {
    mu_coord <- L
    return(list(mu = mu_coord, flag = flag, msg = msg))
  } else if (is.finite(F0) && F0 < alpha) {
    flag <- 1L
    msg  <- "Inactive constraint at 0; returning 0."
    mu_coord <- L
    return(list(mu = mu_coord, flag = flag, msg = msg))
  }
  
  ## F0 > alpha: expand U to bracket the root
  U <- L + U_s
  FU <- Fgamma(U, muA_fixed, muB_fixed)
  while (is.finite(FU) && FU > alpha && U < U_max) {
    U  <- U * U_f
    FU <- Fgamma(U, muA_fixed, muB_fixed)
  }
  if (!is.finite(FU) || FU > alpha) {
    flag <- 1L
    msg  <- "Failed to bracket root; increase U_max or reduce alpha."
    mu_coord <- L
    return(list(mu = mu_coord, flag = flag, msg = msg))
  }
  
  ## Bisection phase
  for (j in 1:MaxIter_b) {
    mid <- L + (U - L) / 2
    if ((U - L) / 2 < delta) break
    Fm <- Fgamma(mid, muA_fixed, muB_fixed)
    if (!is.finite(Fm)) { L <- mid; next }
    if (Fm > alpha) L <- mid else U <- mid
  }
  mu_coord <- L + (U - L) / 2
  list(mu = mu_coord, flag = flag, msg = msg)
}

## -----------------------------
## Algorithm 1: Main Solver for Pi3 Objective
## -----------------------------
## Coordinate Descent to find the optimal vector mu*
solve_mu_for_df_Pi3 <- function(df, Qgrid, init_mu=c(0,0,0)) {
  cache <- make_df_cache(df, Qgrid)
  mu    <- as.numeric(init_mu)
  
  ## Fgamma wrappers matching the subroutine’s signature
  Fwrap_FWER0 <- function(mu1, mu2) {
    function(x, muA, muB) {
      R <- Ri_all_Pi3(c(x, mu1, mu2), cache)
      fwer_values_from_R(R, cache)[["FWER0"]]
    }
  }
  Fwrap_FWER1 <- function(mu0, mu2) {
    function(x, muA, muB) {
      R <- Ri_all_Pi3(c(mu0, x, mu2), cache)
      fwer_values_from_R(R, cache)[["FWER1"]]
    }
  }
  Fwrap_FWER2 <- function(mu0, mu1) {
    function(x, muA, muB) {
      R <- Ri_all_Pi3(c(mu0, mu1, x), cache)
      fwer_values_from_R(R, cache)[["FWER2"]]
    }
  }
  
  for (t in 1:max_outer) {
    old <- mu
    
    ## Update μ0 (Controls FWER0)
    F0fun  <- Fwrap_FWER0(mu[2], mu[3])
    F0_at0 <- F0fun(0, NA, NA)
    mu[1]  <- if (is.finite(F0_at0) && F0_at0 <= alpha) 0 else
      ComputeCoordinateMu(F0fun, muA_fixed = mu[2], muB_fixed = mu[3], alpha = alpha)$mu
    
    ## Update μ1 (Controls FWER1)
    F1fun  <- Fwrap_FWER1(mu[1], mu[3])
    F1_at0 <- F1fun(0, NA, NA)
    mu[2]  <- if (is.finite(F1_at0) && F1_at0 <= alpha) 0 else
      ComputeCoordinateMu(F1fun, muA_fixed = mu[1], muB_fixed = mu[3], alpha = alpha)$mu
    
    ## Update μ2 (Controls FWER2)
    F2fun  <- Fwrap_FWER2(mu[1], mu[2])
    F2_at0 <- F2fun(0, NA, NA)
    mu[3]  <- if (is.finite(F2_at0) && F2_at0 <= alpha) 0 else
      ComputeCoordinateMu(F2fun, muA_fixed = mu[1], muB_fixed = mu[2], alpha = alpha)$mu
    
    if (DBG) {
      Cvals <- fwer_values_from_R(Ri_all_Pi3(mu, cache), cache)
      cat(sprintf("  iter %d: mu=(%.6g, %.6g, %.6g), C=(%.4f, %.4f, %.4f)\n",
                  t, mu[1], mu[2], mu[3], Cvals[1], Cvals[2], Cvals[3]))
    }
    if (sqrt(sum((mu - old)^2)) < outer_tol) break
  }
  list(mu=mu, cache=cache)
}

## -----------------------------
## Classical procedures (Unchanged)
## -----------------------------
bonf_reject    <- function(p, alpha=alpha) as.integer(p <= alpha / length(p))
holm_reject    <- function(p, alpha=alpha) as.integer(p.adjust(p, "holm")     <= alpha)
hoch_reject    <- function(p, alpha=alpha) as.integer(p.adjust(p, "hochberg") <= alpha)
hommel_reject  <- function(p, alpha=alpha) as.integer(p.adjust(p, "hommel")   <= alpha)

## Romano–Wolf stepdown (independence, minP closed form)
romano_wolf_stepdown_indep <- function(p, alpha=alpha) {
  ord <- order(p)              # stepdown in increasing p
  ps  <- p[ord]
  m   <- length(ps)
  padj <- numeric(m)
  for (i in 1:m) {
    mi <- m - i + 1
    padj[i] <- 1 - (1 - ps[i])^mi
  }
  padj <- cummax(padj)           # enforce monotone stepdown
  rej_ord <- as.integer(padj <= alpha)
  rej <- integer(m); rej[ord] <- rej_ord
  rej
}

## Our method: apply Π3-optimized policy to one p-row
our_method_reject_Pi3_t <- function(p_row, mu, df) {
  ord <- order(p_row); ps <- p_row[ord]
  
  ## Compute g(u) values for the observed p-values
  g1 <- g_from_u_t(ps[1], df)
  g2 <- g_from_u_t(ps[2], df)
  g3 <- g_from_u_t(ps[3], df)
  a  <- 2 * g1 * g2 * g3
  
  ## R_i logic is unchanged
  R1 <- a - 6*mu[1]                 - 2*(g2 + g3)*mu[2] - 2*(g2*g3)*mu[3]
  R2 <- a                         - 2*g1*mu[2]        - 2*(g1*g3)*mu[3]
  R3 <- a                                               - 2*(g1*g2)*mu[3]
  
  S <- c(0, R1, R1+R2, R1+R2+R3)
  lstar <- which.max(S) - 1L
  
  rej <- integer(3)
  if (lstar > 0) rej[ord[seq_len(lstar)]] <- 1L
  rej
}

## -----------------------------
## Power estimation under h_3 (all three false)
## -----------------------------
estimate_power_all_t <- function(df, mu, nrep=nrep) {
  # Sample X from the t-distribution
  X  <- matrix(rg_t(nrep * K, df), ncol=K)
  # Convert X to p-values using the N(0,1) null
  P  <- 2 * pnorm(-abs(X))
  
  R_bonf <- rowSums(t(apply(P, 1, bonf_reject,   alpha=alpha)))
  R_holm <- rowSums(t(apply(P, 1, holm_reject,   alpha=alpha)))
  R_hoch <- rowSums(t(apply(P, 1, hoch_reject,   alpha=alpha)))
  R_homm <- rowSums(t(apply(P, 1, hommel_reject, alpha=alpha)))
  R_rw   <- rowSums(t(apply(P, 1, romano_wolf_stepdown_indep, alpha=alpha)))
  R_our  <- rowSums(t(apply(P, 1, our_method_reject_Pi3_t,      mu=mu, df=df)))
  
  list(
    Pi3 = c(
      Bonferroni = mean(R_bonf)/K,
      Holm       = mean(R_holm)/K,
      Hochberg   = mean(R_hoch)/K,
      Hommel     = mean(R_homm)/K,
      RomanoWolf = mean(R_rw)/K,
      OurMethod  = mean(R_our)/K
    ),
    Piany = c(
      Bonferroni = mean(R_bonf > 0),
      Holm       = mean(R_holm > 0),
      Hochberg   = mean(R_hoch > 0),
      Hommel     = mean(R_homm > 0),
      RomanoWolf = mean(R_rw > 0),
      OurMethod  = mean(R_our > 0)
    )
  )
}

## -----------------------------
## MAIN: solve μ (Π3 objective), then simulate & build tables
## -----------------------------
Qgrid <- sample_Q(Ngrid, K=K)

res_list <- list()
mu_warm <- c(0,0,0)

for (df_val in dfs) {
  cat(sprintf("Solving μ (OMT-Π3, t-dist) for df = %.1f ...\n", df_val))
  sol <- solve_mu_for_df_Pi3(df=df_val, Qgrid=Qgrid, init_mu=mu_warm)
  mu_hat <- sol$mu;  mu_warm <- mu_hat
  print(mu_hat)
  
  cat("Simulating power for all procedures (t-dist)...\n")
  est <- estimate_power_all_t(df=df_val, mu=mu_hat, nrep=nrep)
  
  res_list[[as.character(df_val)]] <- list(df=df_val, mu=mu_hat, Pi3=est$Pi3, Piany=est$Piany)
}

## Table A: Π_df,3 (average power) --- policy optimized for Π3
tabA <- do.call(rbind, lapply(res_list, function(z) c(df=z$df, z$Pi3)))
tabA <- as.data.frame(tabA); tabA$df <- as.numeric(tabA$df)
tabA[,-1] <- lapply(tabA[,-1], as.numeric)

## Table B: Π_df,any (probability of ≥1 rejection) --- same policy
tabB <- do.call(rbind, lapply(res_list, function(z) c(df=z$df, z$Piany)))
tabB <- as.data.frame(tabB); tabB$df <- as.numeric(tabB$df)
tabB[,-1] <- lapply(tabB[,-1], as.numeric)

fmt <- function(x) sprintf("%.3f", x)

cat("\n=== Table A: Average power Π_df,3 (policy optimized for Π3; t-dist model) ===\n")
print(data.frame(
  df         = tabA$df,
  Bonferroni = fmt(tabA$Bonferroni),
  Holm       = fmt(tabA$Holm),
  Hochberg   = fmt(tabA$Hochberg),
  Hommel     = fmt(tabA$Hommel),
  RomanoWolf = fmt(tabA$RomanoWolf),
  OurMethod  = fmt(tabA$OurMethod)
), row.names = FALSE)

cat("\n=== Table B: Π_df,any (same Π3-optimized policy; t-dist model) ===\n")
print(data.frame(
  df         = tabB$df,
  Bonferroni = fmt(tabB$Bonferroni),
  Holm       = fmt(tabB$Holm),
  Hochberg   = fmt(tabB$Hochberg),
  Hommel     = fmt(tabB$Hommel),
  RomanoWolf = fmt(tabB$RomanoWolf),
  OurMethod  = fmt(tabB$OurMethod)
), row.names = FALSE)







# --- Start of Plotting Code ---

# 1. Load required libraries
library(ggplot2)
library(tidyr)
library(patchwork)

# 2. Make 'Method' a factor to control order
method_levels <- c("Bonferroni", "Holm", "Hochberg", "Hommel", "RomanoWolf", "OurMethod")
method_labels <- c("Bonferroni", "Holm", "Hochberg", "Hommel", "RomanoWolf", "Algorithm 1")

# 3. Pivot Table A (using the 'tabA' data frame from the simulation)
tabA_long <- pivot_longer(tabA, cols = -df, names_to = "Method", values_to = "Power")
tabA_long$Method <- factor(tabA_long$Method, levels = method_levels, labels = method_labels)

# 4. Pivot Table B (using the 'tabB' data frame from the simulation)
tabB_long <- pivot_longer(tabB, cols = -df, names_to = "Method", values_to = "Power")
tabB_long$Method <- factor(tabB_long$Method, levels = method_levels, labels = method_labels)

# 5. --- MODIFIED Plot 1 (Pi_3) ---
p1 <- ggplot(tabA_long, aes(x = df, y = Power, color = Method, shape = Method, linetype = Method)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 3, alpha = 0.8) +
  labs(title = NULL,
       x = "Degrees of Freedom (df)",
       y = expression(paste("Average Power (", Pi[3], ")"))) +
  scale_x_continuous(breaks = dfs) +  # Shows natural order
  # --- THIS LINE IS CHANGED ---
  scale_y_continuous(limits = c(0, 0.4), breaks = seq(0, 0.4, 0.1), labels = function(x) sprintf("%.1f", x)) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5), 
    legend.title = element_blank(),
    legend.position = c(0.8, 0.6), # Adjusted position for new zoom
    legend.background = element_rect(fill = "white", colour = "grey80")
  )

# 6. --- MODIFIED Plot 2 (Pi_any) ---
p2 <- ggplot(tabB_long, aes(x = df, y = Power, color = Method, shape = Method, linetype = Method)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 3, alpha = 0.8) +
  labs(title = NULL,
       x = "Degrees of Freedom (df)",
       y = expression(paste("Minimal Power (", Pi[any], ")"))) +
  scale_x_continuous(breaks = dfs) + # Shows natural order
  # --- THIS LINE IS CHANGED ---
  scale_y_continuous(limits = c(0, 0.4), breaks = seq(0, 0.4, 0.1), labels = function(x) sprintf("%.1f", x)) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5), sc
    legend.position = "none"
  )

# 7. --- MODIFIED Combination ---
combined_plot <- p1 + p2

# 8. Display and save the final plot
print(combined_plot)
ggsave("path to directory", combined_plot, width = 12, height = 6)
# --- End of Plotting Code ---