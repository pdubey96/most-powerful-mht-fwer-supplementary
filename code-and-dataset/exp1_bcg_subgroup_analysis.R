# ==============================================================================
# Script Name: exp1_bcg_subgroup_analysis.R
# Purpose:     Application of the Optimal Multiple Testing Procedure (Algorithm 1)
#              to the BCG Vaccine Trial Dataset.
#              Reproduces the real-data results (Section 5.1) and Figure 3.
#
# Dataset Details:
#   Name:      'dat.bcg' from the 'metadat' R package.
#   Source:    Colditz, G. A., et al. (1994). Efficacy of BCG vaccine in the 
#              prevention of tuberculosis: Meta-analysis of the published literature. 
#              JAMA, 271(9), 698-702.
#   Content:   Results from 13 clinical trials examining BCG vaccine efficacy.
#
#   Variables Used:
#     - tpos/tneg: TB cases/non-cases in the treated (vaccinated) arm.
#     - cpos/cneg: TB cases/non-cases in the control (non-vaccinated) arm.
#     - alloc:     Method of treatment allocation (random, alternate, systematic).
#     - ablat:     Absolute latitude of the study site (degrees).
#
#   Constructed Subgroups (Splits):
#     The analysis partitions trials into K=3 disjoint subgroups based on:
#     1. Allocation: Random vs Alternate vs Systematic.
#     2. Latitude:   Low vs Mid vs High (tertiles of 'ablat').
#     3. Risk:       Low vs Mid vs High (tertiles of control-arm risk).
#     4. Era:        Early vs Mid vs Late (tertiles of study date).
#
# Methodology:
#   1. Compute fixed-effects estimates (log risk ratio) and z-scores per subgroup.
#   2. Apply Algorithm 1 (Coordinate Descent) to find optimal dual variables (mu).
#   3. Compare decisions with Holm, Hommel, and Closed-Stouffer procedures.
#   4. Generate plots visualizing the "Baseline Risk" finding.
# ==============================================================================

suppressPackageStartupMessages({
  if (!requireNamespace("metafor", quietly = TRUE)) install.packages("metafor")
  if (!requireNamespace("metadat", quietly = TRUE)) install.packages("metadat")
  # parallel is optional; only used if use_parallel <- TRUE
  if (!requireNamespace("parallel", quietly = TRUE)) install.packages("parallel")
})
library(metafor); library(metadat)

# ------------------------------------------------------------------------------
# CONFIGURATION
# ------------------------------------------------------------------------------
alpha            <- 0.05
options(digits = 16)

# Monte Carlo sizes for integration over the ordered simplex Q
# Coarse grid for bracketing; Fine grid for final refinement.
MC_COARSE_SIZE   <- 5e4
MC_FINE_SIZE     <- 1e5

# Solver tolerances and limits for Algorithm 1
TOL              <- 1e-3   # Bisection tolerance
EPS              <- 5e-4   # Outer loop convergence tolerance
TMAX             <- 30     # Max outer iterations

# Beta(a,1) clamp for density estimation (ensures monotonicity, Lemma 2.8)
A_CLAMP_LOW      <- 0.10
A_CLAMP_HIGH     <- 0.90

library(parallel)
# Panel options
RUN_PANEL        <- TRUE                 # Run extra 3-way splits beyond 'alloc'
PANEL_SPLITS     <- c("alloc","lat","risk","era")  # Splits to analyze
GLOBAL_G_PANEL   <- FALSE                 # If TRUE: estimate global 'a' from all data
USE_PARALLEL     <- FALSE                # Parallelize panel splits
N_CORES          <- max(1, parallel::detectCores() - 1)

# File outputs
OUT_SINGLE_TRIPLE_CSV <- "bcg_three_pvals.csv"
OUT_SINGLE_RESULTS_CSV<- "bcg_optionA_results.csv"
OUT_PANEL_RESULTS_CSV <- "bcg_panel_results.csv"
OUT_PANEL_AGG_CSV     <- "bcg_panel_aggregates.csv"

set.seed(1)

# ------------------------------------------------------------------------------
# HELPER FUNCTIONS: Meta-Analysis
# ------------------------------------------------------------------------------
# Fixed-effects meta-analysis for a data subset
get_summary <- function(df) {
  fit <- rma.uni(yi, vi, data = df, method = "FE")
  c(beta = as.numeric(fit$b[1, 1]), se = fit$se)
}

# Calculate z-score and 2-sided p-value
pool_zp <- function(df) {
  s <- get_summary(df)
  z <- s["beta"]/s["se"]
  p <- 2*pnorm(-abs(z))
  list(z = as.numeric(z), p = as.numeric(p))  # list (so [[ "z" ]] works)
}

# ------------------------------------------------------------------------------
# HELPER FUNCTIONS: Density Estimation (g(u))
# ------------------------------------------------------------------------------
# MLE for parameter 'a' in Beta(a,1) model for p-values
est_a_from_p <- function(p, lo = A_CLAMP_LOW, hi = A_CLAMP_HIGH) {
  p <- p[p > 0 & p < 1]
  a_hat <- -length(p) / sum(log(p))
  min(max(a_hat, lo), hi)
}
# Density function closure
g_beta_a1 <- function(a) { force(a); function(u) a * u^(a - 1) }

# ------------------------------------------------------------------------------
# HELPER FUNCTIONS: Monte Carlo & Caching
# ------------------------------------------------------------------------------
# Generate uniform samples on the ordered simplex Q
sample_Q <- function(N) {
  M <- matrix(runif(3*N), ncol = 3)
  t(apply(M, 1, sort))
}

# Pre-compute density values on the grid for fast integration
build_cache <- function(UQ, g_fun) {
  g1 <- g_fun(UQ[,1]); g2 <- g_fun(UQ[,2]); g3 <- g_fun(UQ[,3])
  A  <- 2 * g1 * g2 * g3 # The objective weight function a_i(u)
  list(
    UQ  = UQ, g_fun = g_fun,
    g1  = g1, g2 = g2, g3 = g3,
    A   = A,
    g2p3= g2 + g3,
    g2g3= g2 * g3,
    g1g3= g1 * g3,
    g1g2= g1 * g2
  )
}

# ------------------------------------------------------------------------------
# CORE ALGORITHM: Dual Variables & Decisions
# ------------------------------------------------------------------------------
# Compute R_i functions (see manuscript)
R_all <- function(mu, cache) {
  with(cache, {
    R1 <- A - 6*mu[1] - 2*mu[2]*g2p3 - 2*mu[3]*g2g3
    R2 <- A - 2*mu[2]*g1     - 2*mu[3]*g1g3
    R3 <- A - 2*mu[3]*g1g2
    list(R1=R1,R2=R2,R3=R3)
  })
}

# Determine optimal indicators alpha_i (see manuscript)
alpha_vec <- function(mu, cache) {
  R <- R_all(mu, cache)
  a1 <- (R$R1 > 0) | ((R$R1 + R$R2) > 0) | ((R$R1 + R$R2 + R$R3) > 0)
  a2 <- (R$R2 > 0) | ((R$R2 + R$R3) > 0)
  a3 <- (R$R3 > 0)
  list(a1 = a1, a2 = a2, a3 = a3)
}

# Compute constraint function values F_gamma (see manuscript)
F_maps_vec <- function(mu, cache) {
  a <- alpha_vec(mu, cache)
  F0 <- 6 * mean(a$a1)
  F1 <- 1 - 2 * mean((!a$a2) * cache$g1)                # beta_2 = 1 - alpha_2
  F2 <- 2 * mean(a$a3 * cache$g1 * cache$g2)
  c(F0=F0, F1=F1, F2=F2)
}

# Wrappers for root-finding per coordinate
F0_map <- function(mu0, mu1, mu2, cache) F_maps_vec(c(mu0,mu1,mu2), cache)["F0"]
F1_map <- function(mu1, mu0, mu2, cache) F_maps_vec(c(mu0,mu1,mu2), cache)["F1"]
F2_map <- function(mu2, mu0, mu1, cache) F_maps_vec(c(mu0,mu1,mu2), cache)["F2"]

# ------------------------------------------------------------------------------
# ROOT FINDING: Bracketing & Bisection (Algorithm 2)
# ------------------------------------------------------------------------------
bracket_root <- function(Ffun, alpha, start_L=0, U0=1, grow=2, Umax=1e6, max_expand=40) {
  L <- start_L; fL <- Ffun(L) - alpha
  if (is.na(fL)) stop("F(L) is NA.")
  if (abs(fL) <= 1e-6 || fL < 0) return(c(L=0, U=0))
  U <- U0; fU <- Ffun(U) - alpha; k <- 0L
  while (fU > 0 && U < Umax && k < max_expand) {
    U <- U * grow; fU <- Ffun(U) - alpha; k <- k + 1L
  }
  if (fU > 0) stop(sprintf("Bracket failed: F(0)=%.6g, F(%g)=%.6g (>α).", fL+alpha, U, fU+alpha))
  c(L=L, U=U)
}

bisect_root <- function(Ffun, alpha, L, U, tol=TOL, maxit=80) {
  if (U == 0) return(0)
  for (i in 1:maxit) {
    M  <- 0.5*(L+U); fM <- Ffun(M) - alpha
    if (is.na(fM)) stop("F(M) is NA.")
    if (abs(fM) <= tol || (U-L) < 1e-8) return(M)
    if (fM > 0) L <- M else U <- M
  }
  warning("Bisection maxed; midpoint returned."); 0.5*(L+U)
}

# Coordinate update routine
compute_mu_coord <- function(which, mu_cur, alpha, cache, tol=TOL) {
  if (which == 0) {
    Ffun <- function(x) F0_map(x, mu_cur[2], mu_cur[3], cache)
  } else if (which == 1) {
    Ffun <- function(x) F1_map(x, mu_cur[1], mu_cur[3], cache)
  } else {
    Ffun <- function(x) F2_map(x, mu_cur[1], mu_cur[2], cache)
  }
  br <- bracket_root(Ffun, alpha, start_L=0, U0=1, grow=2, Umax=1e6, max_expand=40)
  bisect_root(Ffun, alpha, br["L"], br["U"], tol=tol, maxit=80)
}

# ------------------------------------------------------------------------------
# MAIN SOLVER: Algorithm 1 (Coordinate Descent)
# ------------------------------------------------------------------------------
solve_mu_optimal <- function(alpha, cache_coarse, cache_fine, eps=EPS, tol=TOL, Tmax=TMAX) {
  mu <- c(0,0,0)
  # Coarse passes (fast convergence)
  for (t in 1:ceiling(Tmax/2)) {
    mu0 <- compute_mu_coord(0, mu, alpha, cache_coarse, tol)
    mu1 <- compute_mu_coord(1, c(mu0, mu[2], mu[3]), alpha, cache_coarse, tol)
    mu2 <- compute_mu_coord(2, c(mu0, mu1, mu[3]), alpha, cache_coarse, tol)
    new <- c(mu0,mu1,mu2)
    if (sqrt(sum((new-mu)^2)) <= eps) { mu <- new; break }
    mu <- new
  }
  # Refinement pass on fine cache
  mu0 <- compute_mu_coord(0, mu, alpha, cache_fine, tol)
  mu1 <- compute_mu_coord(1, c(mu0, mu[2], mu[3]), alpha, cache_fine, tol)
  mu2 <- compute_mu_coord(2, c(mu0, mu1, mu[3]), alpha, cache_fine, tol)
  c(mu0,mu1,mu2)
}

# Apply decisions to a specific sorted p-value vector
opt_decisions_sorted <- function(u_sorted, mu_hat, g_fun) {
  cache1 <- build_cache(matrix(u_sorted, nrow=1), g_fun)
  a <- alpha_vec(mu_hat, cache1)
  as.logical(c(a$a1, a$a1 & a$a2, a$a1 & a$a2 & a$a3))
}

# ------------------------------------------------------------------------------
# BASELINE PROCEDURES
# ------------------------------------------------------------------------------
holm_reject <- function(p, alpha = 0.05) {
  p_ord <- sort(p); K <- length(p)
  crit  <- alpha / (K:1)
  k_fail <- which(p_ord > crit)
  rej_ord <- if (length(k_fail) == 0) rep(TRUE, K) else c(rep(TRUE, min(k_fail) - 1), rep(FALSE, K - min(k_fail) + 1))
  rej <- rep(FALSE, K); rej[order(p)] <- rej_ord
  rej
}
hommel_reject <- function(p, alpha = 0.05) p.adjust(p, method = "hommel") <= alpha

stouffer_local <- function(p_sub, alpha = 0.05) {
  z_ind <- qnorm(1 - p_sub/2) * sign(0.5 - p_sub)
  Zm <- sum(z_ind) / sqrt(length(p_sub))
  abs(Zm) >= qnorm(1 - alpha/2)
}
closed_stouffer_reject <- function(p, alpha = 0.05) {
  stopifnot(length(p) == 3)
  subsets <- list(c(1), c(2), c(3), c(1,2), c(1,3), c(2,3), c(1,2,3))
  rejS <- sapply(subsets, function(S) stouffer_local(p[S], alpha))
  sapply(1:3, function(i) all(rejS[sapply(subsets, function(S) i %in% S)]))
}

# ------------------------------------------------------------------------------
# DATA LOADING
# ------------------------------------------------------------------------------
dat <- metadat::dat.bcg
dat_es <- escalc(measure = "RR", ai = tpos, bi = tneg, ci = cpos, di = cneg, data = dat)

# ------------------------------------------------------------------------------
# SINGLE OUTCOME: Allocation Split
# ------------------------------------------------------------------------------
alloc_fac   <- factor(dat$alloc)
alloc_lvls  <- levels(alloc_fac)
out_alloc   <- lapply(alloc_lvls, function(L) pool_zp(dat_es[ alloc_fac == L, ]))
z_alloc     <- sapply(out_alloc, `[[`, "z"); names(z_alloc) <- alloc_lvls
p_alloc     <- sapply(out_alloc, `[[`, "p"); names(p_alloc) <- alloc_lvls

# Sorted triple for Algorithm 1
ord <- order(p_alloc); inv <- order(ord)
u_sorted <- sort(p_alloc)

# Fit g(u) and build caches
a_hat_single <- est_a_from_p(u_sorted)
g_fun_single <- g_beta_a1(a_hat_single)

UQ_coarse <- sample_Q(MC_COARSE_SIZE)
UQ_fine   <- sample_Q(MC_FINE_SIZE)
cache_coarse <- build_cache(UQ_coarse, g_fun_single)
cache_fine   <- build_cache(UQ_fine,   g_fun_single)

# Solve for mu* and apply decisions
mu_hat_single <- solve_mu_optimal(alpha, cache_coarse, cache_fine, eps=EPS, tol=TOL, Tmax=TMAX)
d_sorted      <- opt_decisions_sorted(u_sorted, mu_hat_single, g_fun_single)
d_opt_single  <- d_sorted[inv]; names(d_opt_single) <- names(p_alloc)

# Apply baselines
rej_holm    <- holm_reject(p_alloc, alpha)
rej_hommel  <- hommel_reject(p_alloc, alpha)
rej_closedS <- closed_stouffer_reject(p_alloc, alpha)

# Residuals (for check)
res_vec_single <- F_maps_vec(mu_hat_single, cache_fine)

res_alloc <- data.frame(
  subgroup = names(p_alloc),
  p_value  = p_alloc,
  reject_Holm = rej_holm,
  reject_Hommel = rej_hommel,
  reject_ClosedStouffer = rej_closedS,
  reject_Optimal_Thm3_8 = d_opt_single
)
cat("\n=== alloc split (single outcome) ===\n")
print(res_alloc, row.names = FALSE)
cat(sprintf("\nmu_hat = (%.6f, %.6f, %.6f),  a_hat = %.4f\n", mu_hat_single[1], mu_hat_single[2], mu_hat_single[3], a_hat_single))
cat(sprintf("constraint residuals: |F0-α|=%.3g, |F1-α|=%.3g, |F2-α|=%.3g\n",
            abs(res_vec_single["F0"]-alpha), abs(res_vec_single["F1"]-alpha), abs(res_vec_single["F2"]-alpha)))

# Save output files
write.csv(
  data.frame(outcome_id = "BCG_alloc", u1 = u_sorted[1], u2 = u_sorted[2], u3 = u_sorted[3]),
  OUT_SINGLE_TRIPLE_CSV, row.names = FALSE
)
write.csv(res_alloc, OUT_SINGLE_RESULTS_CSV, row.names = FALSE)

# ------------------------------------------------------------------------------
# PANEL ANALYSIS: Extra 3-way Splits
# ------------------------------------------------------------------------------
if (RUN_PANEL) {
  # define splits
  splits <- list()
  
  if ("alloc" %in% PANEL_SPLITS) {
    splits$alloc <- factor(dat$alloc)
  }
  if ("lat" %in% PANEL_SPLITS) {
    lat_cut <- cut(dat$ablat, breaks = quantile(dat$ablat, probs = c(0,1/3,2/3,1), na.rm=TRUE),
                   include.lowest = TRUE, labels = c("low_lat","mid_lat","high_lat"))
    splits$lat <- factor(lat_cut)
  }
  if ("risk" %in% PANEL_SPLITS) {
    cer <- with(dat, cpos / (cpos + cneg))
    risk_cut <- cut(cer, breaks = quantile(cer, probs = c(0,1/3,2/3,1), na.rm=TRUE),
                    include.lowest = TRUE, labels = c("low_risk","mid_risk","high_risk"))
    splits$risk <- factor(risk_cut)
  }
  if ("era" %in% PANEL_SPLITS) {
    idx <- seq_len(nrow(dat))
    era_cut <- cut(idx, breaks = quantile(idx, probs = c(0,1/3,2/3,1)),
                   include.lowest = TRUE, labels = c("era1","era2","era3"))
    splits$era <- factor(era_cut)
  }
  
  # Function to run ONE 3-way split
  run_threeway_split <- function(split_fac, split_name,
                                 g_fun_global=NULL, cache_coarse_global=NULL, cache_fine_global=NULL) {
    split_fac <- droplevels(split_fac)
    if (nlevels(split_fac) != 3) return(NULL)
    levs <- levels(split_fac)
    # ensure each level has at least one study
    counts <- table(split_fac)
    if (any(counts == 0)) return(NULL)
    
    out <- lapply(levs, function(L) pool_zp(dat_es[ split_fac == L, ]))
    z <- sapply(out, `[[`, "z"); names(z) <- levs
    p <- sapply(out, `[[`, "p"); names(p) <- levs
    
    ord <- order(p); inv <- order(ord)
    u_sorted <- sort(p)
    
    # choose g: global or local
    if (!is.null(g_fun_global)) {
      g_fun_loc      <- g_fun_global
      cache_coarse   <- cache_coarse_global
      cache_fine     <- cache_fine_global
    } else {
      a_hat_loc      <- est_a_from_p(u_sorted)
      g_fun_loc      <- g_beta_a1(a_hat_loc)
      cache_coarse   <- build_cache(sample_Q(MC_COARSE_SIZE), g_fun_loc)
      cache_fine     <- build_cache(sample_Q(MC_FINE_SIZE),   g_fun_loc)
    }
    
    mu_hat_loc <- solve_mu_optimal(alpha, cache_coarse, cache_fine, eps=EPS, tol=TOL, Tmax=TMAX)
    d_sorted   <- opt_decisions_sorted(u_sorted, mu_hat_loc, g_fun_loc)
    d_opt      <- d_sorted[inv]; names(d_opt) <- names(p)
    
    # baselines
    rej_holm    <- holm_reject(p, alpha)
    rej_hommel  <- hommel_reject(p, alpha)
    rej_closedS <- closed_stouffer_reject(p, alpha)
    
    data.frame(
      outcome_id = split_name,
      level      = names(p),
      p_value    = p,
      reject_Holm = rej_holm,
      reject_Hommel = rej_hommel,
      reject_ClosedStouffer = rej_closedS,
      reject_Optimal_Thm3_8 = d_opt,
      mu0 = mu_hat_loc[1], mu1 = mu_hat_loc[2], mu2 = mu_hat_loc[3],
      a_hat = if (is.null(g_fun_global)) est_a_from_p(u_sorted) else a_hat_single  # report something sensible
    )
  }
  
  # Optionally estimate a single global g from ALL panel triples
  g_fun_global <- NULL; cache_coarse_global <- NULL; cache_fine_global <- NULL
  if (GLOBAL_G_PANEL) {
    # build ALL panel triples first (rough) to estimate a globally
    triples <- list()
    for (nm in names(splits)) {
      fac <- splits[[nm]]
      levs <- levels(fac)
      # skip malformed splits
      if (nlevels(fac) != 3 || any(table(fac) == 0)) next
      tmp <- lapply(levs, function(L) pool_zp(dat_es[ fac == L, ]))
      p   <- sapply(tmp, `[[`, "p")
      triples[[nm]] <- sort(p)
    }
    # include alloc triple too
    triples[["alloc"]] <- sort(p_alloc)
    all_p <- unlist(triples)
    a_hat_global <- est_a_from_p(all_p)
    g_fun_global <- g_beta_a1(a_hat_global)
    cache_coarse_global <- build_cache(sample_Q(MC_COARSE_SIZE), g_fun_global)
    cache_fine_global   <- build_cache(sample_Q(MC_FINE_SIZE),   g_fun_global)
    cat(sprintf("\n[Panel] Using GLOBAL g(u)=Beta(a,1) with a_hat = %.4f\n", a_hat_global))
  }
  
  # run splits
  run_one <- function(nm) {
    run_threeway_split(
      split_fac = splits[[nm]],
      split_name = nm,
      g_fun_global = g_fun_global,
      cache_coarse_global = cache_coarse_global,
      cache_fine_global   = cache_fine_global
    )
  }
  
  res_list <- if (USE_PARALLEL) {
    cl <- parallel::makeCluster(N_CORES)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::parLapply(cl, names(splits), run_one)
  } else {
    lapply(names(splits), run_one)
  }
  panel <- do.call(rbind, res_list)
  panel <- panel[!is.na(panel$outcome_id), , drop=FALSE]
  
  if (nrow(panel) > 0) {
    write.csv(panel, OUT_PANEL_RESULTS_CSV, row.names = FALSE)
    
    # aggregates by outcome_id: avg #rejections, fraction with ≥1 (per-method)
    agg <- do.call(rbind, lapply(split(panel, panel$outcome_id), function(df) {
      sum_by_level <- aggregate(cbind(reject_Holm, reject_Hommel, reject_ClosedStouffer, reject_Optimal_Thm3_8) ~ level, df, sum)
      any_by_level <- aggregate(cbind(reject_Holm, reject_Hommel, reject_ClosedStouffer, reject_Optimal_Thm3_8) ~ level, df, any)
      data.frame(
        outcome_id = df$outcome_id[1],
        avg_num_rej_Holm    = mean(sum_by_level$reject_Holm),
        avg_num_rej_Hommel  = mean(sum_by_level$reject_Hommel),
        avg_num_rej_ClosedS = mean(sum_by_level$reject_ClosedStouffer),
        avg_num_rej_Optimal = mean(sum_by_level$reject_Optimal_Thm3_8),
        frac_any_rej_Holm    = mean(any_by_level$reject_Holm),
        frac_any_rej_Hommel  = mean(any_by_level$reject_Hommel),
        frac_any_rej_ClosedS = mean(any_by_level$reject_ClosedStouffer),
        frac_any_rej_Optimal = mean(any_by_level$reject_Optimal_Thm3_8)
      )
    }))
    write.csv(agg, OUT_PANEL_AGG_CSV, row.names = FALSE)
    
    cat("\n=== Panel summary (head) ===\n")
    print(head(panel, 12), row.names = FALSE)
    cat("\n=== Aggregates ===\n")
    print(agg, row.names = FALSE)
  } else {
    cat("\n[Panel] No valid 3-way splits produced results (check data/levels).\n")
  }
}

cat("\nSaved files:\n",
    sprintf("- %s (alloc triple, sorted)\n", OUT_SINGLE_TRIPLE_CSV),
    sprintf("- %s (alloc split, method comparisons)\n", OUT_SINGLE_RESULTS_CSV),
    if (RUN_PANEL) sprintf("- %s (panel results)\n- %s (panel aggregates)\n", OUT_PANEL_RESULTS_CSV, OUT_PANEL_AGG_CSV) else ""
)

# ------------------------------------------------------------------------------
# POST-PROCESSING: Disagreement Analysis
# ------------------------------------------------------------------------------
panel <- read.csv("bcg_panel_results.csv")

# Identify rows where Optimal differs from baselines
disagree <- subset(panel,
                   reject_Optimal_Thm3_8 != reject_Hommel |
                     reject_Optimal_Thm3_8 != reject_ClosedStouffer |
                     reject_Optimal_Thm3_8 != reject_Holm)
# Display disagreements sorted by outcome and p-value
cat("\n=== Disagreements between Methods ===\n")
print(disagree[order(disagree$outcome_id, disagree$p_value),])

# ------------------------------------------------------------------------------
# DEEP DIVE: Risk Split (Manual verification of g and KKT conditions)
# ------------------------------------------------------------------------------
# Re-construct the 'risk' split manually to inspect variables
risk_fac <- with(dat, {
  cer <- cpos / (cpos + cneg)
  cut(cer, breaks = quantile(cer, probs = c(0,1/3,2/3,1), na.rm=TRUE),
      include.lowest = TRUE, labels = c("low_risk","mid_risk","high_risk"))
})
risk_fac <- factor(risk_fac)

levs <- levels(risk_fac)
get_p <- function(L) {
  s <- get_summary(dat_es[risk_fac == L, ])
  2*pnorm(-abs(s["beta"]/s["se"]))
}
p_risk <- setNames(sapply(levs, get_p), levs)
print(p_risk)

# Solve specifically for this outcome
ord <- order(p_risk); inv <- order(ord)
u_sorted <- sort(p_risk)

a_hat_loc <- est_a_from_p(u_sorted)
g_fun_loc <- g_beta_a1(a_hat_loc)
cache_coarse <- build_cache(sample_Q(MC_COARSE_SIZE), g_fun_loc)
cache_fine   <- build_cache(sample_Q(MC_FINE_SIZE),   g_fun_loc)

mu_hat_loc <- solve_mu_optimal(alpha, cache_coarse, cache_fine, eps=EPS, tol=TOL, Tmax=TMAX)

# Final decisions
d_sorted <- opt_decisions_sorted(u_sorted, mu_hat_loc, g_fun_loc)
d_opt    <- setNames(d_sorted[inv], names(p_risk))

# Baselines for comparison
rej_holm    <- holm_reject(p_risk, alpha)
rej_hommel  <- hommel_reject(p_risk, alpha)
rej_closedS <- closed_stouffer_reject(p_risk, alpha)

# Constraint residuals (F values)
res_vec <- F_maps_vec(mu_hat_loc, cache_fine)

cat("\n=== Risk Split Deep Dive ===\n")
print(cbind(
  p_value = p_risk,
  Holm = rej_holm,
  Hommel = rej_hommel,
  ClosedStouffer = rej_closedS,
  Optimal_Thm3_8 = d_opt
))
cat(sprintf("\n[risk split] a_hat=%.4f, mu_hat=(%.6f, %.6f, %.6f)\n", a_hat_loc, mu_hat_loc[1], mu_hat_loc[2], mu_hat_loc[3]))
cat(sprintf("constraint residuals: |F0-α|=%.3g, |F1-α|=%.3g, |F2-α|=%.3g\n",
            abs(res_vec["F0"]-alpha), abs(res_vec["F1"]-alpha), abs(res_vec["F2"]-alpha)))

# ------------------------------------------------------------------------------
# VERIFICATION: KKT Conditions / Complementary Slackness
# ------------------------------------------------------------------------------
# 1) Print raw F-values
print(res_vec)

# 2) Check binding status
bind_tol <- 5e-3
cat("\nBinding checks:\n")
cat(sprintf("F0~α? %s  (F0=%.4f)\n", abs(res_vec['F0'] - alpha) <= bind_tol, res_vec['F0']))
cat(sprintf("F1≤α?  %s  (F1=%.4f)\n", res_vec['F1'] <= alpha + bind_tol, res_vec['F1']))
cat(sprintf("F2~α? %s  (F2=%.4f)\n", abs(res_vec['F2'] - alpha) <= bind_tol, res_vec['F2']))

# 3) Complementary slackness: products mu * (constraint_slack) should be ~0
kkt <- c(
  mu_hat_loc[1] * (res_vec['F0'] - alpha),   # equality constraint
  mu_hat_loc[2] * (alpha - res_vec['F1']),   # inequality (<=)
  mu_hat_loc[3] * (res_vec['F2'] - alpha)    # equality constraint
)
names(kkt) <- c("mu0*(F0-α)","mu1*(α-F1)","mu2*(F2-α)")
print(kkt)

# ------------------------------------------------------------------------------
# PLOTTING CODE (Produces Manuscript Figures)
# ------------------------------------------------------------------------------
# --- Start of Final Plotting Code ---

# 1. Load required libraries
library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork)
library(scales) # Required for log scale formatting

# 2. Define standard method names and colors for consistency
method_levels <- c("Holm", "Hommel", "ClosedStouffer", "OurMethod")
method_labels <- c("Holm", "Hommel", "Closed-Stouffer", "Algorithm 1")
rejection_colors <- c("Rejected" = "#00B258", "Not Rejected" = "#F04431")

# 3. Filter the 'panel' data frame for the "risk" outcome
risk_data_long <- panel %>%
  filter(outcome_id == "risk") %>%
  select(level, p_value, starts_with("reject_")) %>%
  rename(
    Holm = reject_Holm,
    Hommel = reject_Hommel,
    ClosedStouffer = reject_ClosedStouffer,
    OurMethod = reject_Optimal_Thm3_8
  ) %>%
  pivot_longer(
    cols = c(Holm, Hommel, ClosedStouffer, OurMethod),
    names_to = "Method",
    values_to = "Rejected"
  ) %>%
  mutate(
    Method = factor(Method, levels = method_levels, labels = method_labels),
    level = factor(level, levels = c("low_risk", "mid_risk", "high_risk"),
                   labels = c("Low Risk", "Mid Risk", "High Risk")), # This sets "Low Risk" first
    RejectionStatus = factor(ifelse(Rejected, "Rejected", "Not Rejected"))
  )

# 4. --- Create Facet Labels ---
# Create a data frame for the new facet labels
facet_label_data <- distinct(risk_data_long, level, p_value) %>%
  mutate(
    p_formatted = case_when(
      p_value < 0.001 ~ format(p_value, digits = 3, scientific = TRUE),
      TRUE ~ as.character(round(p_value, 3))
    ),
    facet_label = paste0(level, " (p \u2248 ", p_formatted, ")")
  ) %>%
  select(level, facet_label)

# Join the new facet labels back to the main data
risk_data_long <- risk_data_long %>%
  left_join(facet_label_data, by = "level")

# Re-order the new 'facet_label' factor
risk_data_long$facet_label <- factor(
  risk_data_long$facet_label,
  levels = (facet_label_data %>% arrange(level))$facet_label
)

# 5. Create the new plot
focused_plot <- ggplot(risk_data_long, aes(x = Method, y = p_value)) +
  
  # --- Facet by the level name directly ---
  facet_wrap(~ level) +
  
  # Add points, colored and shaped by rejection status
  geom_point(
    aes(fill = RejectionStatus, shape = RejectionStatus),
    size = 7,      
    stroke = 1.5   
  ) +
  
  # Log scale for y-axis with standard math notation
  scale_y_log10(
    limits = c(1e-30, 2.0), # Ensure points and labels are visible
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  
  # Manual colors and shapes for rejection status
  scale_fill_manual(values = rejection_colors, name = "Decision") +
  scale_shape_manual(values = c("Rejected" = 21, "Not Rejected" = 21), name = "Decision") +
  
  # Remove Title and Subtitle
  labs(
    x = NULL, # Remove x-axis title
    y = "p-value (Log Scale)"
  ) +
  
  theme_bw(base_size = 18) + 
  theme(
    plot.title = element_blank(),  
    plot.subtitle = element_blank(),
    strip.background = element_rect(fill = "grey90", color = "black"),
    strip.text = element_text(size = 16), 
    legend.position = "bottom",
    legend.title = element_text(size = 16, face = "bold"), 
    legend.text = element_text(size = 14), 
    panel.grid.major.x = element_blank(),  
    panel.grid.minor.x = element_blank(),
    axis.title.y = element_text(size = 16), 
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14), 
    axis.text.y = element_text(size = 12) 
  )

# 6. Save and Print
ggsave("path to directory", focused_plot, width = 12, height = 5.5)

print(focused_plot)

# --- End of Focused Plotting Code ---