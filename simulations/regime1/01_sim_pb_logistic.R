## ==============================================================================
## File: 01_sim_pb_logistic.R
## Description: Regime I - Perturbation Bootstrap (PB) for Logistic Regression
## ==============================================================================

suppressPackageStartupMessages({
  library(CVXR)
  library(MASS)
  suppressWarnings(try(library(brglm), silent = TRUE))
  suppressWarnings(try(library(brglm2), silent = TRUE))
})

## ---------------------------------------------------------
## Simulation Configurations
## ---------------------------------------------------------
set.seed(2026)
d_param    <- 10          # Dimension d in {10, 25, 40}
n_sample   <- 50          # Sample size n in {50, 150, 300, 500}
M_mc       <- 500         # Monte Carlo replications (500 or 1000)
B_boot     <- n_sample    # Bootstrap resamples (set B = n)
alpha_nom  <- 0.10        # Nominal significance level (90% CI)

## ---------------------------------------------------------
## Ground Truth & Fixed Design Matrix
## ---------------------------------------------------------
beta_true <- -0.25 + 0.5 * sqrt(1:d_param) * ((-1)^(1:d_param))

# Fixed correlated Gaussian design: Cov(X_j, X_k) = 0.1^|j-k|
mu_zero <- rep(0, d_param)
cov_mat <- outer(1:d_param, 1:d_param, function(i, j) 0.1^(abs(i - j)))
X_raw   <- mvrnorm(n_sample, mu = mu_zero, Sigma = cov_mat)
Z_fixed <- scale(X_raw, center = TRUE, scale = TRUE)

# Fixed success probabilities
prob_true <- 1 / (1 + exp(-as.numeric(Z_fixed %*% beta_true)))

## ---------------------------------------------------------
## Solver: Pure Logistic MLE via CVXR
## ---------------------------------------------------------
solve_logit_mle <- function(y_vec, Z, warm = NULL) {
  n <- nrow(Z); d <- ncol(Z)
  b <- Variable(d)
  if (!is.null(warm)) {
    try({ b$value <- matrix(as.numeric(warm), nrow = d) }, silent = TRUE)
  }
  
  z_b <- Z %*% b
  # Negative log-likelihood: -sum(y * z_b) + sum(log(1 + exp(z_b)))
  nll <- -(1 / n) * sum_entries(matrix(y_vec, n, 1) * z_b) + (1 / n) * sum_entries(logistic(z_b))
  prob <- Problem(Minimize(nll))
  
  for (sv in c("ECOS", "SCS")) {
    res <- tryCatch({
      if (sv == "ECOS") {
        solve(prob, solver = "ECOS", feastol = 1e-8, reltol = 1e-8, abstol = 1e-8, max_iters = 2000, verbose = FALSE)
      } else {
        solve(prob, solver = "SCS", eps = 2e-4, max_iters = 10000, verbose = FALSE)
      }
    }, error = function(e) NULL)
    
    if (!is.null(res) && (res$status %in% c("optimal", "optimal_inaccurate"))) {
      val <- tryCatch(as.numeric(res$getValue(b)), error = function(e) rep(NA_real_, d))
      if (!any(is.na(val))) return(list(ok = TRUE, beta_hat = val))
    }
  }
  list(ok = FALSE, beta_hat = rep(NA_real_, d))
}

## ---------------------------------------------------------
## Output Trackers
## ---------------------------------------------------------
ci_width_mat <- matrix(NA_real_, M_mc, d_param)
cov_ts_mat   <- matrix(NA_integer_, M_mc, d_param)
cov_rs_mat   <- matrix(NA_integer_, M_mc, d_param)
cov_norm_vec <- rep(NA_integer_, M_mc)

cat(sprintf("Running PB Logistic: d = %d, n = %d, M = %d\n", d_param, n_sample, M_mc))

## ---------------------------------------------------------
## Monte Carlo Replication Loop
## ---------------------------------------------------------
for (m in 1:M_mc) {
  y_obs <- rbinom(n_sample, size = 1, prob = prob_true)
  
  warm_coef <- NULL
  if ("brglm" %in% loadedNamespaces()) {
    suppressWarnings({
      fit_br <- try(brglm::brglm.fit(Z_fixed, y_obs, family = binomial(link = "logit")), silent = TRUE)
      if (!inherits(fit_br, "try-error")) warm_coef <- as.numeric(fit_br$coefficients)
    })
  }
  
  fit_mle <- solve_logit_mle(y_obs, Z_fixed, warm = warm_coef)
  if (!fit_mle$ok) fit_mle <- solve_logit_mle(y_obs, Z_fixed, warm = NULL)
  if (!fit_mle$ok) next
  
  beta_hat <- fit_mle$beta_hat
  t_stat_norm <- sqrt(n_sample) * sqrt(sum((beta_hat - beta_true)^2))
  
  prob_hat <- 1 / (1 + exp(-as.numeric(Z_fixed %*% beta_hat)))
  
  beta_boot <- matrix(NA_real_, B_boot, d_param)
  b_idx <- 1
  
  while (b_idx <= B_boot) {
    G_star  <- rexp(n_sample, rate = 1)
    weights <- (y_obs - prob_hat) * (2 - G_star)
    
    v <- Variable(d_param)
    try({ v$value <- matrix(beta_hat, nrow = d_param) }, silent = TRUE)
    
    z_v <- Z_fixed %*% v
    nll_boot  <- -(1 / n_sample) * sum_entries(matrix(y_obs, n_sample, 1) * z_v) + 
                  (1 / n_sample) * sum_entries(logistic(z_v))
    corr_boot <-  (1 / n_sample) * sum_entries(multiply(matrix(weights, n_sample, 1), z_v))
    prob_b    <- Problem(Minimize(nll_boot + corr_boot))
    
    res_b <- tryCatch({
      solve(prob_b, solver = "ECOS", feastol = 1e-8, reltol = 1e-8, abstol = 1e-8, max_iters = 2000, verbose = FALSE)
    }, error = function(e) NULL)
    
    if (is.null(res_b) || !(res_b$status %in% c("optimal", "optimal_inaccurate"))) {
      res_b <- tryCatch({
        solve(prob_b, solver = "SCS", eps = 2e-4, max_iters = 10000, verbose = FALSE)
      }, error = function(e) NULL)
    }
    
    if (!is.null(res_b) && (res_b$status %in% c("optimal", "optimal_inaccurate"))) {
      v_val <- tryCatch(as.numeric(res_b$getValue(v)), error = function(e) rep(NA_real_, d_param))
      if (!any(is.na(v_val))) {
        beta_boot[b_idx, ] <- v_val
        b_idx <- b_idx + 1
      }
    }
  }
  
  tau_boot <- sqrt(n_sample) * (beta_boot - matrix(beta_hat, nrow = B_boot, ncol = d_param, byrow = TRUE))
  tau_norm <- apply(tau_boot, 1, function(r) sqrt(sum(r^2)))
  
  q_low   <- apply(tau_boot, 2, quantile, probs = alpha_nom / 2,     type = 8)
  q_high  <- apply(tau_boot, 2, quantile, probs = 1 - alpha_nom / 2, type = 8)
  q_right <- apply(tau_boot, 2, quantile, probs = alpha_nom,         type = 8)
  q_norm  <- quantile(tau_norm, probs = 1 - alpha_nom, type = 8)
  
  lower_ci    <- beta_hat - (q_high / sqrt(n_sample))
  upper_ci    <- beta_hat - (q_low / sqrt(n_sample))
  right_limit <- beta_hat - (q_right / sqrt(n_sample))
  
  cov_ts_mat[m, ]   <- as.integer(beta_true >= lower_ci & beta_true <= upper_ci)
  cov_rs_mat[m, ]   <- as.integer(beta_true <= right_limit)
  ci_width_mat[m, ] <- (upper_ci - lower_ci)
  cov_norm_vec[m]   <- as.integer(t_stat_norm <= q_norm)
  
  if (m %% 100 == 0) cat(sprintf("  Completed Monte Carlo replication %d / %d\n", m, M_mc))
}

## ---------------------------------------------------------
## Empirical Summaries & MCSE
## ---------------------------------------------------------
p_hat_ts <- colMeans(cov_ts_mat, na.rm = TRUE)
p_hat_rs <- colMeans(cov_rs_mat, na.rm = TRUE)
avg_w    <- colMeans(ci_width_mat, na.rm = TRUE)

mcse_ts_bar <- mean(sqrt(p_hat_ts * (1 - p_hat_ts) / M_mc))
mcse_rs_bar <- mean(sqrt(p_hat_rs * (1 - p_hat_rs) / M_mc))
norm_cov_prob <- mean(cov_norm_vec, na.rm = TRUE)

cat("\n================ Summary Results: PB Logistic ================\n")
cat("Two-Sided Coverage (TS):\n"); print(round(p_hat_ts, 4))
cat("Average CI Widths:\n");       print(round(avg_w, 4))
cat("Right-Sided Coverage (RS):\n");print(round(p_hat_rs, 4))
cat(sprintf("Mean MCSE (TS): %.4f | Mean MCSE (RS): %.4f\n", mcse_ts_bar, mcse_rs_bar))
cat(sprintf("Euclidean Norm Region Coverage: %.4f (MCSE: %.4f)\n", 
            norm_cov_prob, sqrt(norm_cov_prob * (1 - norm_cov_prob) / M_mc)))
