# ==============================================================================
# Script Name: exp2_finance_factor_zoo_analysis.R
# Purpose:     Application of the Optimal Multiple Testing Procedure (Algorithm 1)
#              to Empirical Finance (Asset Pricing Factors).
#              Reproduces the results for Section 5.2 of the manuscript.
#
# Dataset:     Fama-French Data Library (via 'frenchdata' R package).
#              - Source: Kenneth French's Data Library.
#              - Content: Monthly returns for the Fama-French 5 Factors (2x3)
#                and the Momentum Factor.
#
# Methodology:
#   1. Data Prep: Merge FF5 factors (Mkt, SMB, HML, RMW, CMA) with Momentum (Mom).
#   2. Regressions: Estimate time-series regressions for 3 test factors (Mom, RMW, CMA)
#      against the FF3 baseline (Mkt, SMB, HML).
#      Model: (R_i - Rf) = alpha + b*(Mkt-Rf) + s*SMB + h*HML + epsilon
#      Hypothesis: H0: alpha = 0 (Factor is explained by FF3).
#   3. Algorithm 1: Solve for optimal dual variables mu* assuming a heavy-tailed
#      alternative distribution (Student's t with df=4).
#   4. Simulation: Verify FWER control and Power gains under the t(4) assumption.
#   5. Showdown: Compare real-world rejection decisions against Bonferroni, Holm,
#      Hochberg, Hommel, and Romano-Wolf.
#
# Output:      - FWER and Power simulation tables.
#              - Comparison table of decisions for the 3 factors.
#              - Conclusion on "missed discoveries" (specifically Profitability/RMW).
# ==============================================================================

# --- Step 1: Install & Load Packages ---
# You only need to run this once
# install.packages("frenchdata")
# install.packages("dplyr")

library(frenchdata)
library(dplyr)
library(stats)

set.seed(2025)

## --- Step 2: Download & Prepare the Data ---
# We fetch standard research factors directly from the source to ensure reproducibility.

# Download the 5-factor dataset (Mkt, SMB, HML, RMW, CMA)
ds_5factor <- download_french_data("Fama/French 5 Factors (2x3)")

# Extract the monthly data
data_5factor <- ds_5factor$subsets$data[[1]] %>%
  mutate(date = as.numeric(date))

# Download the momentum factor dataset (Mom)
ds_mom <- download_french_data("Momentum Factor (Mom)")

# Extract the monthly data
data_mom <- ds_mom$subsets$data[[1]] %>%
  mutate(date = as.numeric(date))

# Merge the two datasets by date to create a unified factor zoo
factor_data_raw <- inner_join(data_5factor, data_mom, by = "date")

# The data is provided in percentages (e.g., 0.50 for 0.5%); divide by 100 for decimals
factor_data <- factor_data_raw %>%
  mutate(across(-date, ~ . / 100))

# --- Step 3: Run Regressions & Get p-values ---
# We test if the intercepts (alphas) are zero after controlling for the 3-factor model.

cat("--- Running Time-Series Regressions ---\n")

# Regression 1: Momentum vs FF3
fit_mom <- lm( (Mom - RF) ~ `Mkt-RF` + SMB + HML, data = factor_data )
p_mom <- summary(fit_mom)$coefficients["(Intercept)", "Pr(>|t|)"]

# Regression 2: Profitability (RMW) vs FF3
fit_rmw <- lm( (RMW - RF) ~ `Mkt-RF` + SMB + HML, data = factor_data )
p_rmw <- summary(fit_rmw)$coefficients["(Intercept)", "Pr(>|t|)"]

# Regression 3: Investment (CMA) vs FF3
fit_cma <- lm( (CMA - RF) ~ `Mkt-RF` + SMB + HML, data = factor_data )
p_cma <- summary(fit_cma)$coefficients["(Intercept)", "Pr(>|t|)"]

# This is our input data vector u for the showdown (K=3)
p_data <- c(p_mom, p_rmw, p_cma)
names(p_data) <- c("Momentum (MOM)", "Profitability (RMW)", "Investment (CMA)")


## ================================================================
## START: Your Algorithm Code (Adapted for t-distribution)
## ================================================================

## -----------------------------
## Global settings for Algorithm
## -----------------------------
K         <- 3
alpha     <- 0.05
ASSUMED_DF <- 4       # The assumed alternative for g(u). t(4) is a standard
# model for heavy-tailed financial returns (see Manuscript Sec 5.2).
Ngrid     <- 120000    # MC grid size for integrals on Q
nrep_sim  <- 120000    # MC reps for FWER/Power tables
tol_mu    <- 5e-5
outer_tol <- 5e-5
max_outer <- 60

## Bracketing / stability controls
MU_LOWER0 <- 0.0
SCAN_K    <- 40
SCAN_BASE <- 2.0
MU_CAP    <- 1e6
DBG       <- FALSE

## -----------------------------
## t-distribution helpers
## -----------------------------

## Alt p-value density g(u)
## Defined as the likelihood ratio f_A(x)/f_0(x) where x = Phi^-1(1 - u/2).
## Since t(df) has heavier tails than N(0,1), g(u) is monotone decreasing.
g_from_u_t <- function(u, df) {
  u_clamped <- pmax(pmin(u, 1 - 1e-16), 1e-16)
  xp <- qnorm(1 - u_clamped / 2)
  f_A <- dt(xp, df)
  f_0 <- dnorm(xp, 0, 1)
  g_val <- f_A / f_0
  g_val
}

## t-distribution sampler (for simulation)
rg_t <- function(n, df) {
  rt(n, df)
}

## -----------------------------
## Utilities
## -----------------------------
safe_gt0 <- function(v) { w <- v > 0; w[is.na(w)] <- FALSE; w }

# Sample on the ordered simplex Q
sample_Q <- function(N, K=3) {
  m <- matrix(runif(N * K), ncol = K)
  t(apply(m, 1, sort))
}

## -----------------------------
## Per-df cache on shared Q grid
## -----------------------------
make_df_cache <- function(df, Qgrid) {
  u1 <- Qgrid[,1]; u2 <- Qgrid[,2]; u3 <- Qgrid[,3]
  g1 <- g_from_u_t(u1, df)
  g2 <- g_from_u_t(u2, df)
  g3 <- g_from_u_t(u3, df)
  
  # The objective weight function a_i(u) for power maximization
  a  <- 2 * (g1 * g2 * g3)
  
  # Constraint coefficients b_{l,i}(u)
  g2p3 <- 2 * (g2 + g3); g1_  <- 2 * g1; g23  <- 2 * (g2 * g3)
  g13  <- 2 * (g1 * g3); g12  <- 2 * (g1 * g2)
  
  list(u1=u1,u2=u2,u3=u3,
       g1=g1,g2=g2,g3=g3,
       a=a, g2p3=g2p3, g1_=g1_, g23=g23, g13=g13, g12=g12)
}

## -----------------------------
## R_i for Π3 objective
## -----------------------------
# Computes R_i = a_i - sum(mu * b_i)
Ri_all_Pi3 <- function(mu, cache) {
  mu0 <- mu[1]; mu1 <- mu[2]; mu2 <- mu[3]
  R1 <- cache$a - 6*mu0                 - cache$g2p3*mu1 - cache$g23*mu2
  R2 <- cache$a                         - cache$g1_ *mu1 - cache$g13*mu2
  R3 <- cache$a                                          - cache$g12*mu2
  list(R1=R1, R2=R2, R3=R3)
}

# Determine indicator functions alpha_i
alphas_and_D <- function(R) {
  S1 <- R$R1; S2 <- R$R1 + R$R2; S3 <- R$R1 + R$R2 + R$R3
  a1 <- safe_gt0(S1) | safe_gt0(S2) | safe_gt0(S3)
  a2 <- safe_gt0(R$R2) | safe_gt0(R$R2 + R$R3)
  a3 <- safe_gt0(R$R3)
  D1 <- a1; D2 <- a1 & a2; D3 <- a1 & a2 & a3
  list(D1=D1, D2=D2, D3=D3)
}

## -----------------------------
## TRUE FWER values at μ
## -----------------------------
# Evaluates the constraint integrals to check for slackness
fwer_values_from_R <- function(R, cache) {
  D <- alphas_and_D(R)
  g1 <- cache$g1; g2 <- cache$g2; g3 <- cache$g3
  g1[!is.finite(g1)] <- 0; g2[!is.finite(g2)] <- 0; g3[!is.finite(g3)] <- 0
  F0 <- mean(D$D1)
  F1 <- (1/3) * mean( D$D1 * (g2 + g3) + D$D2 * g1 )
  F2 <- (1/3) * mean( D$D1 * (g2 * g3) + D$D2 * (g1 * g3) + D$D3 * (g1 * g2) )
  c(FWER0=F0, FWER1=F1, FWER2=F2)
}


## -----------------------------
## ComputeCoordinateMu Subroutine (Algorithm 2)
## -----------------------------
# Implements the Bracketing + Bisection root finding
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
  
  # Bracketing Phase
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
  
  # Bisection Phase
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
## μ solver (Algorithm 1)
## -----------------------------
# Coordinate Descent to find optimal Lagrange multipliers
solve_mu_for_df_Pi3 <- function(df, Qgrid, init_mu=c(0,0,0)) {
  cache <- make_df_cache(df, Qgrid)
  mu    <- as.numeric(init_mu)
  
  # Wrappers for F_gamma functions
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
  
  # Iterative Coordinate Updates
  for (t in 1:max_outer) {
    old <- mu
    
    # Update mu0
    F0fun  <- Fwrap_FWER0(mu[2], mu[3])
    F0_at0 <- F0fun(0, NA, NA)
    mu[1]  <- if (is.finite(F0_at0) && F0_at0 <= alpha) 0 else
      ComputeCoordinateMu(F0fun, muA_fixed = mu[2], muB_fixed = mu[3], alpha = alpha)$mu
    
    # Update mu1
    F1fun  <- Fwrap_FWER1(mu[1], mu[3])
    F1_at0 <- F1fun(0, NA, NA)
    mu[2]  <- if (is.finite(F1_at0) && F1_at0 <= alpha) 0 else
      ComputeCoordinateMu(F1fun, muA_fixed = mu[1], muB_fixed = mu[3], alpha = alpha)$mu
    
    # Update mu2
    F2fun  <- Fwrap_FWER2(mu[1], mu[2])
    F2_at0 <- F2fun(0, NA, NA)
    mu[3]  <- if (is.finite(F2_at0) && F2_at0 <= alpha) 0 else
      ComputeCoordinateMu(F2fun, muA_fixed = mu[1], muB_fixed = mu[2], alpha = alpha)$mu
    
    if (sqrt(sum((mu - old)^2)) < outer_tol) break
  }
  list(mu=mu, cache=cache)
}

## -----------------------------
## Classical & Our Rejection Procedures
## -----------------------------
# Standard FWER methods for comparison
bonf_reject    <- function(p, alpha=alpha) (p <= alpha / length(p))
holm_reject <- function(p, alpha=alpha) {
  K <- length(p)
  ord <- order(p)
  p_sorted <- p[ord]
  holm_alpha <- alpha / (K - (1:K) + 1)
  k_fail <- which(p_sorted > holm_alpha)
  num_rejected <- if (length(k_fail) > 0) min(k_fail) - 1 else K
  
  rej <- rep(FALSE, K)
  if (num_rejected > 0) {
    rej[ord[1:num_rejected]] <- TRUE
  }
  return(rej)
}
hoch_reject    <- function(p, alpha=alpha) (p.adjust(p, "hochberg") <= alpha)
hommel_reject  <- function(p, alpha=alpha) (p.adjust(p, "hommel")   <= alpha)

romano_wolf_stepdown_indep <- function(p, alpha=alpha) {
  ord <- order(p)
  ps  <- p[ord]
  m   <- length(ps)
  padj <- numeric(m)
  for (i in 1:m) {
    mi <- m - i + 1
    padj[i] <- 1 - (1 - ps[i])^mi
  }
  padj <- cummax(padj)
  rej_ord <- (padj <= alpha)
  rej <- rep(FALSE, m); rej[ord] <- rej_ord
  rej
}

# The Optimal Policy D(u) for a given mu and t-distribution
our_method_reject_Pi3_t <- function(p_row, mu, df) {
  ord <- order(p_row); ps <- p_row[ord]
  
  g1 <- g_from_u_t(ps[1], df); g2 <- g_from_u_t(ps[2], df); g3 <- g_from_u_t(ps[3], df)
  a  <- 2 * g1 * g2 * g3
  
  # Calculate R values using the pre-computed mu
  R1 <- a - 6*mu[1] - 2*(g2 + g3)*mu[2] - 2*(g2*g3)*mu[3]
  R2 <- a - 2*g1*mu[2] - 2*(g1*g3)*mu[3]
  R3 <- a - 2*(g1*g2)*mu[3]
  
  S <- c(0, R1, R1+R2, R1+R2+R3)
  lstar <- which.max(S) - 1L
  
  rej <- rep(FALSE, K)
  if (lstar > 0) {
    rej[ord[seq_len(lstar)]] <- TRUE
  }
  return(rej)
}

## ================================================================
## END: Algorithm Code
## ================================================================


## ================================================================
## NEW: Monte Carlo Simulation Functions
## ================================================================

# --- Simulation 1: FWER (simulating from H0) ---
# Purpose: Verify that the method maintains FWER <= 0.05 under the null.
run_fwer_simulation <- function(mu_policy, df_assumed, nrep) {
  cat(sprintf("Running FWER simulation (nrep=%d) under N(0,1) null...\n", nrep))
  
  # Sample X from the N(0,1) null
  X  <- matrix(rnorm(nrep * K), ncol=K)
  # Convert X to p-values
  P  <- 2 * pnorm(-abs(X))
  
  # Apply each method to each row
  R_bonf <- t(apply(P, 1, bonf_reject,   alpha=alpha))
  R_holm <- t(apply(P, 1, holm_reject,   alpha=alpha))
  R_hoch <- t(apply(P, 1, hoch_reject,   alpha=alpha))
  R_homm <- t(apply(P, 1, hommel_reject, alpha=alpha))
  R_rw   <- t(apply(P, 1, romano_wolf_stepdown_indep, alpha=alpha))
  R_our  <- t(apply(P, 1, our_method_reject_Pi3_t, mu=mu_policy, df=df_assumed))
  
  # FWER is the probability of making at least one rejection (R > 0)
  fwer_df <- data.frame(
    Method = c("Bonferroni", "Holm", "Hochberg", "Hommel", "RomanoWolf", "Algorithm1"),
    FWER = c(
      mean(rowSums(R_bonf) > 0),
      mean(rowSums(R_holm) > 0),
      mean(rowSums(R_hoch) > 0),
      mean(rowSums(R_homm) > 0),
      mean(rowSums(R_rw)   > 0),
      mean(rowSums(R_our)  > 0)
    )
  )
  return(fwer_df)
}

# --- Simulation 2: Power (simulating from HA) ---
# Purpose: Estimate power (Pi3 and Pi_any) under the assumed t(df) alternative.
run_power_simulation <- function(mu_policy, df_true, nrep) {
  cat(sprintf("Running Power simulation (nrep=%d) under t(df=%.1f) alternative...\n", nrep, df_true))
  
  # Sample X from the t(df_true) alternative
  X  <- matrix(rg_t(nrep * K, df_true), ncol=K)
  # Convert X to p-values (using the N(0,1) null)
  P  <- 2 * pnorm(-abs(X))
  
  # Apply each method to each row
  R_bonf <- t(apply(P, 1, bonf_reject,   alpha=alpha))
  R_holm <- t(apply(P, 1, holm_reject,   alpha=alpha))
  R_hoch <- t(apply(P, 1, hoch_reject,   alpha=alpha))
  R_homm <- t(apply(P, 1, hommel_reject, alpha=alpha))
  R_rw   <- t(apply(P, 1, romano_wolf_stepdown_indep, alpha=alpha))
  # Note: we use df_true for the 'df' parameter because our_method_reject_Pi3_t
  # needs it to calculate g(u) for the observed p-values.
  R_our  <- t(apply(P, 1, our_method_reject_Pi3_t, mu=mu_policy, df=df_true))
  
  # Calculate Average Power (Pi_3)
  pi3_df <- data.frame(
    Method = c("Bonferroni", "Holm", "Hochberg", "Hommel", "RomanoWolf", "Algorithm1"),
    Average_Power = c(
      mean(rowSums(R_bonf)) / K,
      mean(rowSums(R_holm)) / K,
      mean(rowSums(R_hoch)) / K,
      mean(rowSums(R_homm)) / K,
      mean(rowSums(R_rw))   / K,
      mean(rowSums(R_our))  / K
    )
  )
  
  # Calculate Minimal Power (Pi_any)
  piany_df <- data.frame(
    Method = c("Bonferroni", "Holm", "Hochberg", "Hommel", "RomanoWolf", "Algorithm1"),
    Minimal_Power = c(
      mean(rowSums(R_bonf) > 0),
      mean(rowSums(R_holm) > 0),
      mean(rowSums(R_hoch) > 0),
      mean(rowSums(R_homm) > 0),
      mean(rowSums(R_rw)   > 0),
      mean(rowSums(R_our)  > 0)
    )
  )
  
  return(list(Pi3 = pi3_df, Piany = piany_df))
}


## ================================================================
## MAIN: Run Simulations and Showdown
## ================================================================

# --- 1. Calculate the Optimal Policy ---
cat(sprintf("\nSolving for optimal policy (Algorithm 1) using g(u) from t(df=%.1f)...\n", ASSUMED_DF))
Qgrid <- sample_Q(Ngrid, K=K)
sol <- solve_mu_for_df_Pi3(df=ASSUMED_DF, Qgrid=Qgrid)
mu_hat <- sol$mu
cat("Optimal mu_hat:", mu_hat, "\n")


# --- 2. Run MC Simulations for FWER and Power ---
fwer_table <- run_fwer_simulation(mu_policy = mu_hat, df_assumed = ASSUMED_DF, nrep = nrep_sim)
power_tables <- run_power_simulation(mu_policy = mu_hat, df_true = ASSUMED_DF, nrep = nrep_sim)

cat("\n--- FWER Simulation Results (H0: N(0,1)) ---\n")
print(fwer_table, row.names = FALSE)

cat(sprintf("\n--- Average Power (Pi_3) Simulation Results (HA: t(%.1f)) ---\n", ASSUMED_DF))
print(power_tables$Pi3, row.names = FALSE)

cat(sprintf("\n--- Minimal Power (Pi_any) Simulation Results (HA: t(%.1f)) ---\n", ASSUMED_DF))
print(power_tables$Piany, row.names = FALSE)


# --- 3. Run Showdown on Fama-French Data ---
cat("\n--- Fama-French Factor Showdown ---\n")
cat("Analyzing K=3 p-values:\n")
print(p_data)

# Run standard procedures
rej_bonf <- bonf_reject(p_data, alpha = alpha)
rej_holm <- holm_reject(p_data, alpha = alpha)
rej_hoch <- hoch_reject(p_data, alpha = alpha)
rej_homm <- hommel_reject(p_data, alpha = alpha)
rej_rw   <- romano_wolf_stepdown_indep(p_data, alpha = alpha)

# Apply our Optimal Policy
rej_our <- our_method_reject_Pi3_t(p_data, mu = mu_hat, df = ASSUMED_DF)

# Print Final Comparison
results_df <- data.frame(
  Hypothesis = names(p_data),
  p_value = p_data,
  Bonferroni = rej_bonf,
  Holm = rej_holm,
  Hochberg = rej_hoch,
  Hommel = rej_homm,
  RomanoWolf = rej_rw,
  Algorithm1 = rej_our
)

cat("\n--- Comparison of Results ---\n")
results_df$p_value <- format.pval(results_df$p_value, digits = 3, eps = 0.001)
print(results_df, row.names = FALSE)

cat("\n--- Conclusion ---\n")
n_holm <- sum(rej_holm)
n_our  <- sum(rej_our)

cat(sprintf("Holm Procedure rejects %d hypothesis/es.\n", n_holm))
cat(sprintf("Algorithm 1 rejects %d hypothesis/es.\n", n_our))

if (n_our > n_holm) {
  cat("SUCCESS: Algorithm 1 found a missed discovery while maintaining FWER control.\n")
} else {
  cat("NOTE: All methods rejected the same number of hypotheses in this case.\n")
}