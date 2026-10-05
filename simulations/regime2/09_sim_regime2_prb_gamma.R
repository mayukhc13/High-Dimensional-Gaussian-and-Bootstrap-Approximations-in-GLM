## ==============================================================================
## File: 09_sim_regime2_prb_gamma.R
## Description: Regime II - Pearson's Residual Bootstrap (PRB) for Sparse Gamma Regression
## Standardized Pearson resamples + l1-penalized least squares optimization
## ==============================================================================

suppressPackageStartupMessages({
  library(MASS)
  library(glmnet)
  library(CVXR)
  suppressWarnings(try(library(h2o), silent = TRUE))
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
shape_gam <- 1.0

## ---------------------------------------------------------
## Initialize H2O
## ---------------------------------------------------------
if ("h2o" %in% loadedNamespaces()) {
  tryCatch({
    h2o.init(nthreads = -1, max_mem_size = "4g", verbose = FALSE)
    h2o.no_progress()
  }, error = function(e) message("H2O local cluster connection notice."))
}

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
Z_fixed <- scale(X_raw, center = TRUE, scale = FALSE)

eta_true   <- as.numeric(Z_fixed %*% beta_true)
scale_true <- exp(eta_true) / shape_gam

## ---------------------------------------------------------
## Trackers
## ---------------------------------------------------------
ci_width_mat <- matrix(NA_real_, M_mc, d0_active)
cov_ts_mat   <- matrix(NA_integer_, M_mc, d0_active)
cov_rs_mat   <- matrix(NA_integer_, M_mc, d0_active)
cov_norm_vec <- rep(NA_integer_, M_mc)

cat(sprintf("Running Regime II PRB Gamma: (d, d0) = (%d, %d), n = %d, M = %d\n", 
            d_param, d0_active, n_sample, M_mc))

## ---------------------------------------------------------
## Monte Carlo Replication Loop
## ---------------------------------------------------------
for (m in 1:M_mc) {
  y_obs <- rgamma(n_sample, shape = shape_gam, scale = scale_true)
  
  # Step 1: Lambda via 10-fold CV
  lambda_opt <- NULL
  if (h2o.is_client_connected()) {
    df_gamma <- as.h2o(data.frame(y = y_obs, Z_fixed))
    cv_h2o <- tryCatch({
      h2o.glm(y = "y", x = colnames(df_gamma)[-1], training_frame = df_gamma,
              family = "gamma", link = "log", alpha = 1, nfolds = 10,
              standardize = FALSE, intercept = FALSE)
    }, error = function(e) NULL)
    if (!is.null(cv_h2o)) lambda_opt <- h2o.getLambdaBest(cv_h2o)
  }
  
  if (is.null(lambda_opt) || is.na(lambda_opt)) {
    cv_glmnet_approx <- cv.glmnet(Z_fixed, log(pmax(y_obs, 1e-4)), family = "gaussian", 
                                  alpha = 1, nfolds = 10, intercept = FALSE)
    lambda_opt <- cv_glmnet_approx$lambda.min
  }
  
  # Step 2: Base Lasso via CVXR
  v_base    <- Variable(d_param)
  eta_base  <- Z_fixed %*% v_base
  obj_lasso <- (1 / n_sample) * sum_entries(shape_gam * (eta_base + matrix(y_obs, n_sample, 1) * exp(-eta_base))) + 
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
  
  # Threshold screening
  beta_scr <- beta_lasso
  beta_scr[abs(beta_scr) <= (lambda_opt / n_sample)] <- 0
  beta_active_hat <- beta_scr[1:d0_active]
  
  stat_norm <- sqrt(n_sample) * sqrt(sum((beta_active_hat - beta_active_true)^2))
  
  # Step 3: Standardized Pearson Residuals & Matrix G_bar
  # Under Gamma log link:
  # G_bar_i = sqrt(alpha) * x_i
  eta_hat  <- as.numeric(Z_fixed %*% beta_lasso)
  mu_hat   <- exp(eta_hat)
  denom_sd <- mu_hat / sqrt(shape_gam)
  
  e_sharp  <- (y_obs - mu_hat) / denom_sd
  e_center <- e_sharp - mean(e_sharp)
  
  # G_bar matrix: n x d
  G_bar <- sqrt(shape_gam) * Z_fixed
  G_beta_hat <- as.numeric(G_bar %*% beta_lasso)
  
  # Step 4: PRB Resampling via Penalized Least Squares
  beta_boot <- matrix(NA_real_, B_boot, d0_active)
  
  for (b_idx in 1:B_boot) {
    e_star <- sample(e_center, size = n_sample, replace = TRUE)
    y_synth <- G_beta_hat + e_star
    
    fit_prb <- glmnet(G_bar, y_synth, family = "gaussian", intercept = FALSE, 
                      alpha = 1, lambda = lambda_opt / n_sample, standardize = FALSE)
    b_val <- as.numeric(fit_prb$beta)
    
    b_val[abs(b_val) <= (lambda_opt / n_sample)] <- 0
    beta_boot[b_idx, ] <- b_val[1:d0_active]
  }
  
  # Step 5: Quantiles & Coverage Tracking
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

cat("\n================ Results: Regime II PRB Gamma ================\n")
cat("Two-Sided Coverage (TS):\n"); print(round(p_hat_ts, 3))
cat("Average CI Widths:\n");       print(round(avg_w, 2))
cat("Right-Sided Coverage (RS):\n");print(round(p_hat_rs, 3))
cat(sprintf("Mean TS: %.3f | Mean RS: %.3f | Mean MCSE: (%.3f)\n", 
            mean(p_hat_ts), mean(p_hat_rs), mcse_ts_bar))
cat(sprintf("Euclidean Norm Region Coverage: %.3f (MCSE: %.3f)\n", 
            norm_cov_prob, sqrt(norm_cov_prob * (1 - norm_cov_prob) / M_mc)))
