#' Compute NRI and IDI with bootstrap validation
#' @param new_model_result List with lin_pred and bas_cumhazdf
#' @param clin_vars Vector of clinical predictor names
#' @param biomarker_vars Vector of biomarker predictor names
#' @param data Data frame with time, status, predictors, and INDEX column
#' @param eval_time Time point for evaluation
#' @param event_vars Event variable names
#' @param recalibrate Apply calibration slope and intercept to new model predictions
#' @param n_bootstrap Number of bootstrap samples
#' @return Bootstrap results
compute_nri_idi_bootstrap <- function(new_model_result, clin_vars, biomarker_vars, data, eval_time, 
                                      event_vars = c("time", "status"), 
                                      recalibrate = list(slope = 1, intercept = 0), 
                                      n_bootstrap = 1000) {
  # Check and load dependencies
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package 'survival' is required but not installed. Please install it.")
  }
  if (!requireNamespace("survIDINRI", quietly = TRUE)) {
    stop("Package 'survIDINRI' is required but not installed. Please install it with: install.packages('survIDINRI')")
  }
  library(survival)
  library(survIDINRI)
  
  # Standardize covariates
  data_scaled <- data
  all_vars <- c(clin_vars, biomarker_vars)
  for (var in all_vars) {
    data_scaled[[var]] <- scale(data_scaled[[var]], center = TRUE, scale = TRUE)
  }
  
  # Run single dataset computation
  full_result <- compute_nri_idi_single(
    new_model_result = new_model_result,
    clin_vars = clin_vars,
    biomarker_vars = biomarker_vars,
    data = data_scaled,
    eval_time = eval_time,
    event_vars = event_vars,
    recalibrate = recalibrate
  )
  if (is.null(full_result)) {
    stop("Single dataset computation failed")
  }
  
  # Bootstrap loop
  boot_results <- list(
    brier_new = numeric(n_bootstrap),
    brier_ref = numeric(n_bootstrap),
    nri = numeric(n_bootstrap),
    idi = numeric(n_bootstrap)
  )
  warning_count <- 0
  
  for (i in 1:n_bootstrap) {
    set.seed(i)
    boot_idx <- sample(1:nrow(data_scaled), nrow(data_scaled), replace = TRUE)
    boot_data <- data_scaled[boot_idx, ]
    
    boot_result <- compute_nri_idi_single(
      new_model_result = new_model_result,
      clin_vars = clin_vars,
      biomarker_vars = biomarker_vars,
      data = boot_data,
      eval_time = eval_time,
      event_vars = event_vars,
      recalibrate = recalibrate
    )
    
    if (!is.null(boot_result)) {
      boot_results$brier_new[i] <- boot_result$brier_new
      boot_results$brier_ref[i] <- boot_result$brier_ref
      if (!is.null(boot_result$idi_result) && 
          !is.null(boot_result$idi_result$m2.est) && 
          length(boot_result$idi_result$m2.est) >= 1 &&
          !is.null(boot_result$idi_result$m1.est) && 
          length(boot_result$idi_result$m1.est) >= 1) {
        boot_results$nri[i] <- boot_result$idi_result$m2.est[1]  # Continuous NRI
        boot_results$idi[i] <- boot_result$idi_result$m1.est[1]  # IDI
      } else {
        warning_count <- warning_count + 1
        boot_results$nri[i] <- NA
        boot_results$idi[i] <- NA
      }
    } else {
      warning_count <- warning_count + 1
      boot_results$brier_new[i] <- NA
      boot_results$brier_ref[i] <- NA
      boot_results$nri[i] <- NA
      boot_results$idi[i] <- NA
    }
  }
  
  # Summarize warnings
  if (warning_count > 0) {
    cat("Total warnings during bootstrap: ", warning_count, "\n")
  }
  
  # Summarize bootstrap results
  cat("\nBootstrap Internal Validation Results (n =", n_bootstrap, "):\n")
  cat("Brier Score (New Model): Mean =", mean(boot_results$brier_new, na.rm = TRUE), 
      ", 95% CI =", quantile(boot_results$brier_new, c(0.025, 0.975), na.rm = TRUE), "\n")
  cat("Brier Score (Reference Model): Mean =", mean(boot_results$brier_ref, na.rm = TRUE), 
      ", 95% CI =", quantile(boot_results$brier_ref, c(0.025, 0.975), na.rm = TRUE), "\n")
  cat("Continuous NRI: Mean =", mean(boot_results$nri, na.rm = TRUE), 
      ", 95% CI =", quantile(boot_results$nri, c(0.025, 0.975), na.rm = TRUE), "\n")
  cat("IDI: Mean =", mean(boot_results$idi, na.rm = TRUE), 
      ", 95% CI =", quantile(boot_results$idi, c(0.025, 0.975), na.rm = TRUE), "\n")
  
  return(list(
    full_result = full_result,
    boot_results = boot_results
  ))
}

# Example
bootstrap_result <- compute_nri_idi_bootstrap(
  new_model_result = cox_glmnet_result,
  clin_vars = clin_vars3,
  biomarker_vars = top_biomarkers,
  data = npxdata_all,
  eval_time = eval_time,
  event_vars = c("time", "status"),
  recalibrate = list(slope = 2.444079, intercept = 1.756214),
  n_bootstrap = 100
)
