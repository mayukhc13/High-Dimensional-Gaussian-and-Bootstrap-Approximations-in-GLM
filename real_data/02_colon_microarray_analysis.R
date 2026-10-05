## ==============================================================================
## File: 02_colon_microarray_analysis.R
## Description: Real Data Application - Colon Microarray Classification (Logistic Lasso)
## Alon et al. colon dataset (n = 62, d = 2000 genes).
## Performs 10-fold CV, screening, and constructs 90% CIs via PB and PRB.
## ==============================================================================

suppressPackageStartupMessages({
  library(glmnet)
  library(CVXR)
  suppressWarnings(try(library(plsgenomics), silent = TRUE))
})

set.seed(2026)
alpha_nom <- 0.10  # 90% Confidence Intervals

## ---------------------------------------------------------
## Step 1: Load or Replicate Colon Microarray Dataset
## ---------------------------------------------------------
data_loaded <- FALSE

if ("plsgenomics" %in% loadedNamespaces()) {
  tryCatch({
    data(Colon, package = "plsgenomics")
    X_raw <- Colon$X
    y_raw <- Colon$Y
    # Format outcome: 1 = Tumor (positive), 0 = Normal
    y_vec <- ifelse(y_raw == 1, 1, 0)
    data_loaded <- TRUE
    cat("Successfully loaded Colon dataset from 'plsgenomics'.\n")
  }, error = function(e) NULL)
}

if (!data_loaded) {
  cat("Package 'plsgenomics' not detected. Generating structurally identical microarray setup (n=62, d=2000)...\n")
  n_samp <- 62
  p_feat <- 2000
  X_raw  <- matrix(rnorm(n_samp * p_feat), nrow = n_samp, ncol = p_feat)
  colnames(X_raw) <- paste0("Gene_", 1:p_feat)
  
  # Sparse signals on first 5 genes
  signal <- 1.2 * X_raw[, 1] - 1.0 * X_raw[, 2] + 0.8 * X_raw[, 3] - 0.7 * X_raw[, 4] + 0.9 * X_raw[, 5]
  p_syn  <- 1 / (1 + exp(-signal))
  y_vec  <- rbinom(n_samp, size = 1, prob = p_syn)
}

n_obs  <- nrow(X_raw)
d_gene <- ncol(X_raw)
B_boot <- 500

# Standardize gene expressions
Z_design <- scale(X_raw, center = TRUE, scale = TRUE)
gene_names <- colnames(Z_design)
if (is.null(gene_names)) gene_names <- paste0("Gene_", 1:d_gene)

cat(sprintf("Colon Microarray setup: n = %d tissue samples, d = %d genes\n", n_obs, d_gene))

## ---------------------------------------------------------
## Step 2: 10-Fold CV Lasso & Screening
## ---------------------------------------------------------
cv_fit <- cv.glmnet(Z_design, y_vec, family = "binomial", alpha = 1, 
                    nfolds = 10, type.measure = "deviance", intercept = TRUE)
lambda_opt <- cv_fit$lambda.min

fit_lasso  <- glmnet(Z_design, y_vec, family = "binomial", alpha = 1, 
                     lambda = lambda_opt, intercept = TRUE)

beta_lasso <- as.numeric(fit_lasso$beta)
intercept_lasso <- as.numeric(fit_lasso$a0)

# Screening threshold tau_n = lambda / n
active_idx <- which(abs(beta_lasso) > (lambda_opt / n_obs))
d0_active  <- length(active_idx)

cat(sprintf("Optimal lambda = %.5f | Selected %d active genes\n", lambda_opt, d0_active))

if (d0_active == 0) {
  active_idx <- order(abs(beta_lasso), decreasing = TRUE)[1:5]
  d0_active  <- length(active_idx)
  cat(sprintf("Note: Retaining top %d genes for inferential demonstration.\n", d0_active))
}

beta_act_hat <- beta_lasso[active_idx]
active_names <- gene_names[active_idx]

eta_hat  <- as.numeric(intercept_lasso + Z_design %*% beta_lasso)
prob_hat <- 1 / (1 + exp(-eta_hat))

## ---------------------------------------------------------
## Step 3: Perturbation Bootstrap (PB) on Active Support
## ---------------------------------------------------------
cat("Running Perturbation Bootstrap (PB)...\n")
Z_active <- Z_design[, active_idx, drop = FALSE]
beta_pb  <- matrix(NA_real_, B_boot, d0_active)
b_idx    <- 1
retries  <- 0

while (b_idx <= B_boot && retries < 50) {
  G_star  <- rexp(n_obs, rate = 1)
  weights <- (y_vec - prob_hat) * (2 - G_star)
  
  v_b   <- Variable(d0_active)
  eta_b <- Z_active %*% v_b
  
  obj_b <- (1 / n_obs) * sum_entries(logistic(eta_b) - matrix(y_vec, n_obs, 1) * eta_b) + 
           (1 / n_obs) * sum_entries(multiply(matrix(weights, n_obs, 1), eta_b)) + 
           (lambda_opt / n_obs) * p_norm(v_b, 1)
  
  res_b <- tryCatch({
    solve(Problem(Minimize(obj_b)), solver = "ECOS", feastol = 1e-6, abstol = 1e-6, verbose = FALSE)
  }, error = function(e) NULL)
  
  if (!is.null(res_b) && (res_b$status %in% c("optimal", "optimal_inaccurate"))) {
    val <- as.numeric(res_b$getValue(v_b))
    val[abs(val) <= (lambda_opt / n_obs)] <- 0
    beta_pb[b_idx, ] <- val
    b_idx <- b_idx + 1
  } else {
    retries <- retries + 1
  }
}

## ---------------------------------------------------------
## Step 4: Pearson's Residual Bootstrap (PRB)
## ---------------------------------------------------------
cat("Running Pearson's Residual Bootstrap (PRB)...\n")
v_diag   <- pmax(prob_hat * (1 - prob_hat), 1e-5)
e_sharp  <- (y_vec - prob_hat) / sqrt(v_diag)
e_center <- e_sharp - mean(e_sharp)

G_bar      <- sqrt(v_diag) * Z_active
G_beta_hat <- as.numeric(G_bar %*% beta_act_hat)

beta_prb <- matrix(NA_real_, B_boot, d0_active)
for (b in 1:B_boot) {
  e_star  <- sample(e_center, size = n_obs, replace = TRUE)
  y_synth <- G_beta_hat + e_star
  
  fit_b <- glmnet(G_bar, y_synth, family = "gaussian", intercept = FALSE, 
                  alpha = 1, lambda = lambda_opt / n_obs, standardize = FALSE)
  val <- as.numeric(fit_b$beta)
  val[abs(val) <= (lambda_opt / n_obs)] <- 0
  beta_prb[b, ] <- val
}

## ---------------------------------------------------------
## Step 5: Construct 90% Confidence Intervals
## ---------------------------------------------------------
calc_ci <- function(beta_est, boot_mat, n, alpha) {
  tau <- sqrt(n) * (boot_mat - matrix(beta_est, nrow = nrow(boot_mat), ncol = length(beta_est), byrow = TRUE))
  q_lo <- apply(tau, 2, quantile, probs = alpha / 2,     type = 8, na.rm = TRUE)
  q_hi <- apply(tau, 2, quantile, probs = 1 - alpha / 2, type = 8, na.rm = TRUE)
  lower <- beta_est - (q_hi / sqrt(n))
  upper <- beta_est - (q_lo / sqrt(n))
  data.frame(Gene = active_names, Estimate = beta_est, Lower = lower, Upper = upper, Width = upper - lower)
}

ci_pb_df  <- calc_ci(beta_act_hat, beta_pb,  n_obs, alpha_nom)
ci_prb_df <- calc_ci(beta_act_hat, beta_prb, n_obs, alpha_nom)

cat("\n================ Colon Microarray: 90% CI (PB) ================\n")
print(round(ci_pb_df[, -1], 4), row.names = ci_pb_df$Gene)

cat("\n================ Colon Microarray: 90% CI (PRB) ===============\n")
print(round(ci_prb_df[, -1], 4), row.names = ci_prb_df$Gene)
