## ==============================================================================
## File: 05_sim_regime2_gaussian_failure_logistic.R
## Description: Regime II - Finite-Sample Failure of Gaussian Approximation (Logistic)
## Demonstrates severe undercoverage of classical Wald CIs based on active Fisher info
## ==============================================================================

suppressPackageStartupMessages({
  library(MASS)
  library(glmnet)
})

## ---------------------------------------------------------
## Simulation Configurations
## Setups: (d, d0) in {(100,6), (200,8), (300,10), (500,12)}
## ---------------------------------------------------------
set.seed(2026)
d_param   <- 100         # Ambient dimension d
d0_active <- 6           # True sparsity d0
n_sample  <- 100         # Sample size n in {100, 250, 350, 500}
M_mc      <- 500         # Monte Carlo replications (500 or 1000)
alpha_nom <- 0.10        # Nominal 90% Wald CIs (z_0.95 = 1.644853)
z_crit    <- qnorm(1 - alpha_nom / 2)

## ---------------------------------------------------------
## Ground Truth & Fixed Design Matrix
## ---------------------------------------------------------
beta_true <- numeric(d_param)
for (j in 1:d0_active) {
  beta_true[j] <- -0.25 + 0.5 * sqrt(j) * ((-1)^j)
}
beta_active_true <- beta_true[1:d0_active]

# Fixed correlated Gaussian design: Cov(X_j, X_k) = 0.1^|j-k|
cov_mat <- outer(1:d_param, 1:d_param, function(i, j) 0.1^(abs(i - j)))
X_raw   <- mvrnorm(n_sample, mu = rep(0, d_param), Sigma = cov_mat)
Z_fixed <- scale(X_raw, center = TRUE, scale = TRUE)
Z_11    <- Z_fixed[, 1:d0_active, drop = FALSE]

# True probabilities
prob_true <- 1 / (1 + exp(-as.numeric(Z_fixed %*% beta_true)))

## ---------------------------------------------------------
## Trackers
## ---------------------------------------------------------
cov_gauss_mat <- matrix(NA_integer_, M_mc, d0_active)

cat(sprintf("Running Gaussian Breakdown: (d, d0) = (%d, %d), n = %d, M = %d\n", 
            d_param, d0_active, n_sample, M_mc))

## ---------------------------------------------------------
## Monte Carlo Replication Loop
## ---------------------------------------------------------
for (m in 1:M_mc) {
  y_obs <- rbinom(n_sample, size = 1, prob = prob_true)
  
  # 10-fold CV Lasso via cv.glmnet
  cv_fit <- cv.glmnet(Z_fixed, y_obs, family = "binomial", intercept = FALSE, 
                      alpha = 1, nfolds = 10, type.measure = "deviance")
  lambda_opt <- cv_fit$lambda.min
  
  # Extract Lasso estimates (no intercept)
  fit_lasso <- glmnet(Z_fixed, y_obs, family = "binomial", intercept = FALSE, 
                      alpha = 1, lambda = lambda_opt)
  beta_lasso <- as.numeric(fit_lasso$beta)
  
  # Apply conservative screening threshold tau_n = lambda / n
  beta_scr <- beta_lasso
  beta_scr[abs(beta_scr) <= (lambda_opt / n_sample)] <- 0
  beta_active_est <- beta_scr[1:d0_active]
  
  # Fitted probabilities under Lasso estimator
  prob_hat <- 1 / (1 + exp(-as.numeric(Z_fixed %*% beta_lasso)))
  w_diag   <- pmax(prob_hat * (1 - prob_hat), 1e-5)
  
  # Empirical Fisher Information matrix S_n = (1/n) X' W X
  # Sub-matrix S_{n,11} corresponding to active coordinates {1, ..., d0}
  S_n11 <- crossprod(Z_11, w_diag * Z_11) / n_sample
  
  # Invert active Fisher Information
  S_n11_inv <- tryCatch({
    solve(S_n11)
  }, error = function(e) {
    ginv(S_n11)
  })
  
  # Coordinate-wise standard error: sqrt((S_{n,11}^{-1})_{jj} / n)
  se_gauss <- sqrt(pmax(diag(S_n11_inv), 0) / n_sample)
  
  lower_wald <- beta_active_est - z_crit * se_gauss
  upper_wald <- beta_active_est + z_crit * se_gauss
  
  cov_gauss_mat[m, ] <- as.integer(beta_active_true >= lower_wald & beta_active_true <= upper_wald)
  
  if (m %% 100 == 0) cat(sprintf("  Completed Monte Carlo replication %d / %d\n", m, M_mc))
}

## ---------------------------------------------------------
## Empirical Summaries & MCSE
## ---------------------------------------------------------
p_hat_gauss <- colMeans(cov_gauss_mat, na.rm = TRUE)
mcse_gauss  <- sqrt(p_hat_gauss * (1 - p_hat_gauss) / M_mc)

cat("\n================ Results: Gaussian Breakdown (Logistic) ================\n")
for (j in 1:d0_active) {
  cat(sprintf("beta_%-2d: Coverage = %.3f (MCSE = %.3f)\n", j, p_hat_gauss[j], mcse_gauss[j]))
}
cat(sprintf("Average Active Coverage: %.3f (Mean MCSE: %.3f)\n", 
            mean(p_hat_gauss), mean(mcse_gauss)))
