## ==============================================================================
## File: 06_sim_regime2_pb_logistic.R
## Description: Regime II - Perturbation Bootstrap (PB) for Sparse Logistic Regression
## Evaluates coordinate-wise coverage and Euclidean norm confidence regions
## ==============================================================================

suppressPackageStartupMessages({
  library(MASS)
  library(glmnet)
  library(CVXR)
})

## ---------------------------------------------------------
## Simulation Configurations
## Setups: (d, d0) in {(100,6), (200,8), (300,10), (500,12)}
## ---------------------------------------------------------
set.seed(2026)
d_param   <- 100
d0_active <- 6
n_sample  <- 100
M_mc      <- 500
B_boot    <- n_sample
alpha_nom <- 0.10

## ---------------------------------------------------------
## Ground Truth & Fixed Design Matrix
## ---------------------------------------------------------
beta_true <- numeric(d_param)
for (j in 1:d0_active) {
  beta_true[j] <- -0.25 + 0.5 * sqrt(j) * ((-1)^j)
}
beta_active_true <- beta_true[1:d0_active]

cov_mat <- outer(1:d_param, 1:d_param, function(i, j) 0.1^(abs(i - j)))
X_raw   <- mvrnorm(n_sample, mu = rep(0, d_param), Sigma = cov_mat)
Z_fixed <- scale(X_raw, center = TRUE, scale = TRUE)

prob_true <- 1 / (1 + exp(-as.numeric(Z_fixed %*% beta_true)))

## ---------------------------------------------------------
## Trackers
## ---------------------------------------------------------
ci_width_mat <- matrix(NA_real_, M_mc, d0_active)
cov_ts_mat   <- matrix(NA_integer_, M_mc, d0_active)
cov_rs_mat   <- matrix(NA_integer_, M_mc, d0_active)
cov_norm_vec <- rep(NA_integer_, M_mc)

cat(sprintf("Running Regime II PB Logistic: (d, d0) = (%d, %d), n = %d, M = %d\n", 
            d_param, d0_active, n_sample, M_mc))

## ---------------------------------------------------------
## Monte Carlo Replication Loop
## ---------------------------------------------------------
for (m in 1:M_mc) {
  y_obs <- rbinom(n_sample, size = 1, prob = prob_true)
  
  # Step 1: Compute lambda_cv via 10-fold cross validation
  cv_fit <- cv.glmnet(Z_fixed, y_obs, family = "binomial", intercept = FALSE, 
                      alpha = 1, nfolds = 10, type.measure = "deviance")
  lambda_opt <- cv_fit$lambda.min
  
  # Step 2: Solve base Lasso estimator via CVXR
  v_base    <- Variable(d_param)
  eta_base  <- Z_fixed %*% v_base
  obj_lasso <- (1 / n_sample) * sum_entries(logistic(eta_base) - matrix(y_obs, n_sample, 1) * eta_base) + 
               (lambda_opt / n_sample) * p_norm(v_base, 1)
  prob_lasso <- Problem(Minimize(obj_lasso))
  
  res_lasso <- tryCatch({
    solve(prob_lasso, solver = "ECOS", feastol = 1e-7, reltol = 1e-7, abstol = 1e-7, verbose = FALSE)
  }, error = function(e) NULL)
  
  if (is.null(res_lasso) || !(res_lasso$status %in% c("optimal", "optimal_inaccurate"))) {
    res_lasso <- tryCatch({
      solve(prob_lasso, solver = "SCS", eps = 2e-4, max_iters = 10000, verbose = FALSE)
    }, error = function(e) NULL)
  }
  
  if (is.null(res_lasso) || !(res_lasso$status %in% c("optimal", "optimal_inaccurate"))) next
  beta_lasso <- as.numeric(res_lasso$getValue(v_base))
  
  # Finite-sample conservative screening threshold tau_n = lambda / n
  beta_scr <- beta_lasso
  beta_scr[abs(beta_scr) <= (lambda_opt / n_sample)] <- 0
  beta_active_hat <- beta_scr[1:d0_active]
  
  # Statistic for norm-based coverage: || T_n^(1) || = sqrt(n) * || beta_hat^(1) - beta^(1) ||
  stat_norm <- sqrt(n_sample) * sqrt(sum((beta_active_hat - beta_active_true)^2))
  
  # Fitted probability vector at beta_lasso
  eta_hat  <- as.numeric(Z_fixed %*% beta_lasso)
  prob_hat <- 1 / (1 + exp(-eta_hat))
  
  # Step 3: PB Resampling
  beta_boot <- matrix(NA_real_, B_boot, d0_active)
  b_idx <- 1
  retries <- 0
  
  while (b_idx <= B_boot && retries < 30) {
    G_star  <- rexp(n_sample, rate = 1)
    weights <- (y_obs - prob_hat) * (2 - G_star)
    
    v_b <- Variable(d_param)
    try({ v_b$value <- matrix(beta_lasso, nrow = d_param) }, silent = TRUE)
    
    eta_b    <- Z_fixed %*% v_b
    obj_boot <- (1 / n_sample) * sum_entries(logistic(eta_b) - matrix(y_obs, n_sample, 1) * eta_b) + 
                (1 / n_sample) * sum_entries(multiply(matrix(weights, n_sample, 1), eta_b)) + 
                (lambda_opt / n_sample) * p_norm(v_b, 1)
    prob_b   <- Problem(Minimize(obj_boot))
    
    res_b <- tryCatch({
      solve(prob_b, solver = "ECOS", feastol = 1e-7, reltol = 1e-7, abstol = 1e-7, verbose = FALSE)
    }, error = function(e) NULL)
    
    if (is.null(res_b) || !(res_b$status %in% c("optimal", "optimal_inaccurate"))) {
      res_b <- tryCatch({
        solve(prob_b, solver = "SCS", eps = 2e-4, max_iters = 10000, verbose = FALSE)
      }, error = function(e) NULL)
    }
    
    if (!is.null(res_b) && (res_b$status %in% c("optimal", "optimal_inaccurate"))) {
      b_val <- tryCatch(as.numeric(res_b$getValue(v_b)), error = function(e) rep(NA_real_, d_param))
      if (!any(is.na(b_val))) {
        # Screen bootstrap estimates identically
        b_val[abs(b_val) <= (lambda_opt / n_sample)] <- 0
        beta_boot[b_idx, ] <- b_val[1:d0_active]
        b_idx <- b_idx + 1
      } else {
        retries <- retries + 1
      }
    } else {
      retries <- retries + 1
    }
  }
  
  if (b_idx <= B_boot) next
  
  # Step 4: Quantile Calculation on Active Coordinates
  tau_boot <- sqrt(n_sample) * (beta_boot - matrix(beta_active_hat, nrow = B_boot, ncol = d0_active, byrow = TRUE))
  tau_norm <- apply(tau_boot, 1, function(r) sqrt(sum(r^2)))
  
  q_low   <- apply(tau_boot, 2, quantile, probs = alpha_nom / 2,     type = 8, na.rm = TRUE)
  q_high  <- apply(tau_boot, 2, quantile, probs = 1 - alpha_nom / 2, type = 8, na.rm = TRUE)
  q_right <- apply(tau_boot, 2, quantile, probs = alpha_nom,         type = 8, na.rm = TRUE)
  q_norm  <- quantile(tau_norm, probs = 1 - alpha_nom, type = 8, na.rm = TRUE)
  
  lower_ci    <- beta_active_hat - (q_high / sqrt(n_sample))
  upper_ci    <- beta_active_hat - (q_low / sqrt(n_sample))
  right_limit <- beta_active_hat - (q_right / sqrt(n_sample))
  
  cov_ts_mat[m, ]   <- as.integer(beta_active_true >= lower_ci & beta_active_true <= upper_ci)
  cov_rs_mat[m, ]   <- as.integer(beta_active_true <= right_limit)
  ci_width_mat[m, ] <- (upper_ci - lower_ci)
  cov_norm_vec[m]   <- as.integer(stat_norm <= q_norm)
  
  if (m %% 50 == 0) cat(sprintf("  Completed Monte Carlo replication %d / %d\n", m, M_mc))
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

cat("\n================ Results: Regime II PB Logistic ================\n")
cat("Two-Sided Coverage (TS):\n"); print(round(p_hat_ts, 3))
cat("Average CI Widths:\n");       print(round(avg_w, 2))
cat("Right-Sided Coverage (RS):\n");print(round(p_hat_rs, 3))
cat(sprintf("Mean TS: %.3f | Mean RS: %.3f | Mean MCSE: (%.3f)\n", 
            mean(p_hat_ts), mean(p_hat_rs), mcse_ts_bar))
cat(sprintf("Euclidean Norm Region Coverage: %.3f (MCSE: %.3f)\n", 
            norm_cov_prob, sqrt(norm_cov_prob * (1 - norm_cov_prob) / M_mc)))
