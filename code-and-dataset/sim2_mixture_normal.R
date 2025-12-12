## ================================================================
## 2-Component Mixture-Normal setting
## K=3, independent
## H0: X ~ N(0, 1)
## HA: X ~ 0.5 * N(theta, 1) + 0.5 * N(-theta, 1)
## p-values via 2-sided test: p = 2 * pnorm(-abs(X))
## Our method: optimize Π3 (average power) with strong FWER (l=0,1,2)
## Compare: Bonferroni, Holm, Hochberg, Hommel, Romano–Wolf (indep minP)
## Output: Table A (Π3) and Table B (Π_any), policy optimized for Π3
## ================================================================

set.seed(2025)

## -----------------------------
## Global settings
## -----------------------------
K         <- 3
alpha     <- 0.05
# --- THIS LINE IS CHANGED ---
thetas    <- seq(0, -4, by = -0.5)
nrep      <- 120000    # Monte Carlo reps for power
Ngrid     <- 120000    # MC grid size for integrals on Q
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
## Mixture-normal helpers
## -----------------------------
Phi     <- pnorm
PhiInv  <- qnorm

## Alt p-value density g(u)
## p = 2 * pnorm(-abs(X))
## x_p = pnorm(1 - p/2)
## g(p) = LR(x_p) = exp(-theta^2/2) * cosh(x_p * theta)
g_from_u_mixture <- function(u, theta) {
  # Clamp u to avoid Inf at u=0
  u_clamped <- pmax(u, 1e-16)
  xp <- PhiInv(1 - u_clamped / 2)
  
  # Use abs(theta) since mixture is symmetric in theta and -theta
  th_abs <- abs(theta)
  g_val <- exp(-th_abs^2 / 2) * cosh(xp * th_abs)
  
  # When theta=0, g_val should be 1
  g_val[theta == 0] <- 1.0
  g_val
}

## Mixture normal sampler: X ~ 0.5*N(theta,1) + 0.5*N(-theta,1)
rg_mixture <- function(n, theta) {
  # Flip n coins
  coin_flips <- runif(n) < 0.5
  
  # Generate from N(theta, 1)
  sample1 <- rnorm(n, mean = theta, sd = 1)
  
  # Generate from N(-theta, 1)
  sample2 <- rnorm(n, mean = -theta, sd = 1)
  
  # Combine
  ifelse(coin_flips, sample1, sample2)
}


## -----------------------------
## Utilities
## -----------------------------
safe_gt0 <- function(v) { w <- v > 0; w[is.na(w)] <- FALSE; w }

## Monte Carlo grid on Q (order-statistics region)
sample_Q <- function(N, K=3) {
  m <- matrix(runif(N * K), ncol = K)
  t(apply(m, 1, sort))
}

## -----------------------------
## Per-theta cache on shared Q grid
## -----------------------------
make_theta_cache <- function(theta, Qgrid) {
  u1 <- Qgrid[,1]; u2 <- Qgrid[,2]; u3 <- Qgrid[,3]
  g1 <- g_from_u_mixture(u1, theta)
  g2 <- g_from_u_mixture(u2, theta)
  g3 <- g_from_u_mixture(u3, theta)
  
  ## Π3 objective: a_i(u) = 2*g1*g2*g3 for i=1,2,3
  a  <- 2 * (g1 * g2 * g3)
  
  ## Constraint pieces b_{l,i}(u) (from your Lemma; unchanged form)
  g2p3 <- 2 * (g2 + g3)    # for FWER1 with D1
  g1_  <- 2 * g1          # for FWER1 with D2
  g23  <- 2 * (g2 * g3)    # for FWER2 with D1
  g13  <- 2 * (g1 * g3)    # for FWER2 with D2
  g12  <- 2 * (g1 * g2)    # for FWER2 with D3
  
  list(u1=u1,u2=u2,u3=u3,
       g1=g1,g2=g2,g3=g3,
       a=a, g2p3=g2p3, g1_=g1_, g23=g23, g13=g13, g12=g12)
}

## -----------------------------
## R_i for Π3 objective
## -----------------------------
Ri_all_Pi3 <- function(mu, cache) {
  mu0 <- mu[1]; mu1 <- mu[2]; mu2 <- mu[3]
  R1 <- cache$a - 6*mu0                 - cache$g2p3*mu1 - cache$g23*mu2
  R2 <- cache$a                         - cache$g1_ *mu1 - cache$g13*mu2
  R3 <- cache$a                                          - cache$g12*mu2
  list(R1=R1, R2=R2, R3=R3)
}

alphas_and_D <- function(R) {
  S1 <- R$R1
  S2 <- R$R1 + R$R2
  S3 <- R$R1 + R$R2 + R$R3
  a1 <- safe_gt0(S1) | safe_gt0(S2) | safe_gt0(S3)  # α1
  a2 <- safe_gt0(R$R2) | safe_gt0(R$R2 + R$R3)      # α2
  a3 <- safe_gt0(R$R3)                             # α3
  D1 <- a1
  D2 <- a1 & a2
  D3 <- a1 & a2 & a3
  list(D1=D1, D2=D2, D3=D3)
}

## -----------------------------
## TRUE FWER values at μ (unchanged linear forms)
## -----------------------------
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
## ComputeCoordinateMu Subroutine (Unchanged)
## -----------------------------
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
  
  ## F0 > alpha: expand U
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
  
  ## Bisection
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
## μ solver (Algorithm 1) for Π3 objective (Unchanged)
## -----------------------------
solve_mu_for_theta_Pi3 <- function(theta, Qgrid, init_mu=c(0,0,0)) {
  cache <- make_theta_cache(theta, Qgrid)
  mu    <- as.numeric(init_mu)
  
  ## Fgamma wrappers matching the subroutine’s (x, muA, muB) signature
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
    
    ## μ0 (FWER0)
    F0fun  <- Fwrap_FWER0(mu[2], mu[3])
    F0_at0 <- F0fun(0, NA, NA)
    mu[1]  <- if (is.finite(F0_at0) && F0_at0 <= alpha) 0 else
      ComputeCoordinateMu(F0fun, muA_fixed = mu[2], muB_fixed = mu[3], alpha = alpha)$mu
    
    ## μ1 (FWER1)
    F1fun  <- Fwrap_FWER1(mu[1], mu[3])
    F1_at0 <- F1fun(0, NA, NA)
    mu[2]  <- if (is.finite(F1_at0) && F1_at0 <= alpha) 0 else
      ComputeCoordinateMu(F1fun, muA_fixed = mu[1], muB_fixed = mu[3], alpha = alpha)$mu
    
    ## μ2 (FWER2)
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
our_method_reject_Pi3 <- function(p_row, mu, theta) {
  ord <- order(p_row); ps <- p_row[ord]
  
  ## Compute g(u) values for the observed p-values
  g1 <- g_from_u_mixture(ps[1], theta)
  g2 <- g_from_u_mixture(ps[2], theta)
  g3 <- g_from_u_mixture(ps[3], theta)
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
estimate_power_all <- function(theta, mu, nrep=nrep) {
  # Sample X from the mixture
  X  <- matrix(rg_mixture(nrep * K, theta), ncol=K)
  # Convert X to p-values
  P  <- 2 * pnorm(-abs(X))
  
  R_bonf <- rowSums(t(apply(P, 1, bonf_reject,   alpha=alpha)))
  R_holm <- rowSums(t(apply(P, 1, holm_reject,   alpha=alpha)))
  R_hoch <- rowSums(t(apply(P, 1, hoch_reject,   alpha=alpha)))
  R_homm <- rowSums(t(apply(P, 1, hommel_reject, alpha=alpha)))
  R_rw   <- rowSums(t(apply(P, 1, romano_wolf_stepdown_indep, alpha=alpha)))
  R_our  <- rowSums(t(apply(P, 1, our_method_reject_Pi3,      mu=mu, theta=theta)))
  
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

for (th in thetas) {
  cat(sprintf("Solving μ (OMT-Π3, mixture model) for theta = %.2f ...\n", th))
  sol <- solve_mu_for_theta_Pi3(theta=th, Qgrid=Qgrid, init_mu=mu_warm)
  mu_hat <- sol$mu;  mu_warm <- mu_hat
  print(mu_hat)
  
  cat("Simulating power for all procedures (mixture model)...\n")
  est <- estimate_power_all(theta=th, mu=mu_hat, nrep=nrep)
  
  res_list[[as.character(th)]] <- list(theta=th, mu=mu_hat, Pi3=est$Pi3, Piany=est$Piany)
}

## Table A: Π_theta,3 (average power) --- policy optimized for Π3
tabA <- do.call(rbind, lapply(res_list, function(z) c(theta=z$theta, z$Pi3)))
tabA <- as.data.frame(tabA); tabA$theta <- as.numeric(tabA$theta)
tabA[,-1] <- lapply(tabA[,-1], as.numeric)

## Table B: Π_theta,any (probability of ≥1 rejection) --- same policy
tabB <- do.call(rbind, lapply(res_list, function(z) c(theta=z$theta, z$Piany)))
tabB <- as.data.frame(tabB); tabB$theta <- as.numeric(tabB$theta)
tabB[,-1] <- lapply(tabB[,-1], as.numeric)

fmt <- function(x) sprintf("%.3f", x)

cat("\n=== Table A: Average power Π_theta,3 (policy optimized for Π3; mixture model) ===\n")
print(data.frame(
  theta      = tabA$theta,
  Bonferroni = fmt(tabA$Bonferroni),
  Holm       = fmt(tabA$Holm),
  Hochberg   = fmt(tabA$Hochberg),
  Hommel     = fmt(tabA$Hommel),
  RomanoWolf = fmt(tabA$RomanoWolf),
  OurMethod  = fmt(tabA$OurMethod)
), row.names = FALSE)

cat("\n=== Table B: Π_theta,any (same Π3-optimized policy; mixture model) ===\n")
print(data.frame(
  theta      = tabB$theta,
  Bonferroni = fmt(tabB$Bonferroni),
  Holm       = fmt(tabB$Holm),
  Hochberg   = fmt(tabB$Hochberg),
  Hommel     = fmt(tabB$Hommel),
  RomanoWolf = fmt(tabB$RomanoWolf),
  OurMethod  = fmt(tabB$OurMethod)
), row.names = FALSE)







# --- Start of Plotting Code ---

# # 1. Load required libraries
# library(ggplot2)
# library(tidyr)
# library(patchwork)
# 
# # 2. Make 'Method' a factor to control order
# method_levels <- c("Bonferroni", "Holm", "Hochberg", "Hommel", "RomanoWolf", "OurMethod")
# # --- THIS LINE IS ADDED ---
# method_labels <- c("Bonferroni", "Holm", "Hochberg", "Hommel", "RomanoWolf", "Algorithm 1")
# 
# # 3. Pivot Table A (using the 'tabA' data frame from the simulation)
# tabA_long <- pivot_longer(tabA, cols = -theta, names_to = "Method", values_to = "Power")
# # --- THIS LINE IS CHANGED (added labels) ---
# tabA_long$Method <- factor(tabA_long$Method, levels = method_levels, labels = method_labels)
# 
# # 4. Pivot Table B (using the 'tabB' data frame from the simulation)
# tabB_long <- pivot_longer(tabB, cols = -theta, names_to = "Method", values_to = "Power")
# # --- THIS LINE IS CHANGED (added labels) ---
# tabB_long$Method <- factor(tabB_long$Method, levels = method_levels, labels = method_labels)
# 
# # 5. --- MODIFIED Plot 1 (Pi_3) ---
# # We will add the legend *inside* this plot
# p1 <- ggplot(tabA_long, aes(x = theta, y = Power, color = Method, shape = Method, linetype = Method)) +
#   geom_line(linewidth = 0.8) +
#   geom_point(size = 3, alpha = 0.8) +
#   labs(title = NULL,
#        x = expression(paste("Signal Strength (", theta, ")")),
#        y = expression(paste("Average Power (", Pi[3], ")"))) +
#   # --- THIS LINE IS CHANGED (added labels) ---
#   scale_x_reverse(labels = function(x) sprintf("%.1f", x)) + 
#   # --- THIS LINE IS CHANGED (added labels) ---
#   scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), labels = function(x) sprintf("%.1f", x)) +
#   theme_bw(base_size = 14) +
#   theme(
#     plot.title = element_text(hjust = 0.5), 
#     legend.title = element_blank(),
#     legend.position = c(0.22, 0.75), # Place legend inside: 22% from left, 75% from bottom
#     legend.background = element_rect(fill = "white", colour = "grey80") # Add a box
#   )
# 
# # 6. --- MODIFIED Plot 2 (Pi_any) ---
# # We will *remove* the legend from this plot
# p2 <- ggplot(tabB_long, aes(x = theta, y = Power, color = Method, shape = Method, linetype = Method)) +
#   geom_line(linewidth = 0.8) +
#   geom_point(size = 3, alpha = 0.8) +
#   labs(title = NULL,
#        x = expression(paste("Signal Strength (", theta, ")")),
#        y = expression(paste("Minimal Power (", Pi[any], ")"))) +
#   # --- THIS LINE IS CHANGED (added labels) ---
#   scale_x_reverse(labels = function(x) sprintf("%.1f", x)) +
#   # --- THIS LINE IS CHANGED (added labels) ---
#   scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), labels = function(x) sprintf("%.1f", x)) +
#   theme_bw(base_size = 14) +
#   theme(
#     plot.title = element_text(hjust = 0.5), 
#     legend.position = "none" # Remove legend
#   )


# 1. Load required libraries
library(ggplot2)
library(tidyr)
library(patchwork)

# 2. Make 'Method' a factor to control order
method_levels <- c("Bonferroni", "Holm", "Hochberg", "Hommel", "RomanoWolf", "OurMethod")
# --- THIS LINE IS ADDED ---
method_labels <- c("Bonferroni", "Holm", "Hochberg", "Hommel", "RomanoWolf", "Algorithm 1")

# 3. Pivot Table A (using the 'tabA' data frame from the simulation)
tabA_long <- pivot_longer(tabA, cols = -theta, names_to = "Method", values_to = "Power")
# --- THIS LINE IS CHANGED (added labels) ---
tabA_long$Method <- factor(tabA_long$Method, levels = method_levels, labels = method_labels)

# 4. Pivot Table B (using the 'tabB' data frame from the simulation)
tabB_long <- pivot_longer(tabB, cols = -theta, names_to = "Method", values_to = "Power")
# --- THIS LINE IS CHANGED (added labels) ---
tabB_long$Method <- factor(tabB_long$Method, levels = method_levels, labels = method_labels)

# 5. --- MODIFIED Plot 1 (Pi_3) ---
# We will add the legend *inside* this plot
p1 <- ggplot(tabA_long, aes(x = theta, y = Power, color = Method, shape = Method, linetype = Method)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 3, alpha = 0.8) +
  labs(title = NULL,
       x = expression(paste("Signal Strength (", theta, ")")),
       y = expression(paste("Average Power (", Pi[3], ")"))) +
  # --- THIS LINE IS CHANGED (added labels) ---
  scale_x_reverse(labels = function(x) sprintf("%.1f", x)) + 
  # --- THIS LINE IS CHANGED (added labels) ---
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), labels = function(x) sprintf("%.1f", x)) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5), 
    legend.title = element_blank(),
    legend.position = c(0.22, 0.75), # Place legend inside: 22% from left, 75% from bottom
    legend.background = element_rect(fill = "white", colour = "grey80"), # Add a box
    axis.title = element_text(size = 28),
    axis.text  = element_text(size = 24),
    legend.text = element_text(size = 18)
  )

# 6. --- MODIFIED Plot 2 (Pi_any) ---
# We will *remove* the legend from this plot
p2 <- ggplot(tabB_long, aes(x = theta, y = Power, color = Method, shape = Method, linetype = Method)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 3, alpha = 0.8) +
  labs(title = NULL,
       x = expression(paste("Signal Strength (", theta, ")")),
       y = expression(paste("Minimal Power (", Pi[any], ")"))) +
  # --- THIS LINE IS CHANGED (added labels) ---
  scale_x_reverse(labels = function(x) sprintf("%.1f", x)) +
  # --- THIS LINE IS CHANGED (added labels) ---
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), labels = function(x) sprintf("%.1f", x)) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5), 
    legend.position = "none", # Remove legend
    axis.title = element_text(size = 28),
    axis.text  = element_text(size = 24)
  )


# 7. Combine and Save 
p1 <- p1 + theme(plot.margin = margin(5.5, 40, 5.5, 5.5))
p2 <- p2 + theme(plot.margin = margin(5.5, 5.5, 5.5, 40))
combined_plot <- p1 + p2

# 8. Display and save the final plot
print(combined_plot)

# --- CHANGING FILENAME ---
ggsave("path to directory/filename.pdf", combined_plot, width = 13.5, height = 6)
# --- End of Plotting Code ---
