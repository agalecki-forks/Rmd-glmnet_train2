#' Compute NRI and IDI for survival model comparison with internal validation
#' @param new_model_result List with lin_pred and bas_cumhazdf
#' @param clin_vars Vector of clinical predictor names
#' @param biomarker_vars Vector of biomarker predictor names
#' @param data Data frame with time, status, predictors, and INDEX column
#' @param eval_time Time point for evaluation
#' @param event_vars Event variable names
#' @param risk_thresholds Ignored (categorical NRI not supported by IDI.INF)
#' @param recalibrate Apply calibration slope and intercept to new model predictions
#' @param n_bootstrap Number of bootstrap samples
#' @return NRI and IDI results, diagnostics
compute_nri_idi_internal <- function(new_model_result, clin_vars, biomarker_vars, data, eval_time, 
                                    event_vars = c("time", "status"), risk_thresholds = NULL, 
                                    recalibrate = list(slope = 1, intercept = 0), 
                                    n_bootstrap = 1000) {
  # Check and load dependencies
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package 'survival' is required but not installed. Please install it.")
  }
  if (!requireNamespace("survIDINRI", quietly = TRUE)) {
    stop("Package 'survIDINRI' is required but not installed. Please install it with: install.packages('survIDINRI')")
  }
  if (!exists("IDI.INF", envir = asNamespace("survIDINRI"))) {
    stop("Function 'IDI.INF' not found in survIDINRI package. Please ensure the correct package version is installed.")
  }
  library(survival)
  library(survIDINRI)
  
  # Validate inputs
  if (!all(event_vars %in% colnames(data))) {
    stop("Event variables not found in data: ", paste(event_vars, collapse = ", "))
  }
  if (!all(clin_vars %in% colnames(data))) {
    stop("Clinical predictors not found in data: ", paste(clin_vars, collapse = ", "))
  }
  if (!all(biomarker_vars %in% colnames(data))) {
    stop("Biomarker predictors not found in data: ", paste(biomarker_vars, collapse = ", "))
  }
  if (!"INDEX" %in% colnames(data)) {
    stop("INDEX column not found in data")
  }
  time <- data[[event_vars[1]]]
  status <- data[[event_vars[2]]]
  
  # Check for missing values in predictors
  all_vars <- c(clin_vars, biomarker_vars)
  if (any(is.na(data[, all_vars]))) {
    stop("Missing values detected in predictors: ", 
         paste(all_vars[colSums(is.na(data[, all_vars])) > 0], collapse = ", "))
  }
  
  # Check for non-numeric predictors
  non_numeric <- all_vars[!sapply(data[, all_vars], is.numeric)]
  if (length(non_numeric) > 0) {
    stop("Non-numeric predictors detected: ", paste(non_numeric, collapse = ", "))
  }
  
  # Check lin_pred structure
  if (is.null(names(new_model_result$lin_pred)) || length(names(new_model_result$lin_pred)) == 0) {
    warning("names(lin_pred) is NULL or empty. Assuming lin_pred order matches data rows.")
    if (length(new_model_result$lin_pred) != nrow(data)) {
      stop("lin_pred length (", length(new_model_result$lin_pred), 
           ") does not match data rows (", nrow(data), ")")
    }
    names(new_model_result$lin_pred) <- as.character(data$INDEX)
  }
  
  # Create formula for reference model (for Brier score)
  ref_formula <- as.formula(paste("Surv(", event_vars[1], ",", event_vars[2], ") ~", 
                                  paste(clin_vars, collapse = " + ")))
  
  # Function to compute predictions and NRI/IDI for a single dataset
  compute_nri_idi_single <- function(data_subset, new_model_result, ref_formula, eval_time) {
    # Align lin_pred using INDEX column
    lin_pred_new <- new_model_result$lin_pred[match(as.character(data_subset$INDEX), 
                                                    as.character(names(new_model_result$lin_pred)))]
    if (any(is.na(lin_pred_new))) {
      warning("Missing lin_pred values for ", sum(is.na(lin_pred_new)), 
              " observations. Check INDEX alignment.")
      cat("First few INDEX values in data_subset:\n")
      print(head(data_subset$INDEX))
      cat("First few names in lin_pred:\n")
      print(head(names(new_model_result$lin_pred)))
      cat("Non-matching INDEX values:\n")
      print(head(data_subset$INDEX[!data_subset$INDEX %in% names(new_model_result$lin_pred)]))
    }
    
    # Predict survival probabilities for new model
    bas_cumhazdf <- new_model_result$bas_cumhazdf
    cumhaz_at_t <- approx(bas_cumhazdf$time, bas_cumhazdf$hazard, xout = eval_time, 
                          method = "linear", rule = 2)$y
    surv_prob_new <- exp(-cumhaz_at_t * exp(lin_pred_new))
    
    # Recalibrate new model predictions
    if (!is.null(recalibrate$slope) && !is.null(recalibrate$intercept)) {
      logit_pred <- log(pmax(surv_prob_new, 1e-10) / (1 - pmax(surv_prob_new, 1e-10)))
      logit_adj <- recalibrate$intercept + recalibrate$slope * logit_pred
      surv_prob_new <- exp(logit_adj) / (1 + exp(logit_adj))
    }
    risk_new <- 1 - surv_prob_new
    
    # Fit reference model for Brier score
    ref_model <- tryCatch({
      environment(ref_formula) <- list2env(list(data_subset = data_subset), 
                                           parent = .GlobalEnv)
      coxph(ref_formula, data = data_subset, x = TRUE)
    }, error = function(e) {
      warning("Reference model fitting failed: ", e$message)
      return(NULL)
    })
    if (is.null(ref_model)) {
      return(NULL)
    }
    
    # Predict survival probabilities for reference model
    surv_fit_ref <- survfit(ref_model, newdata = data_subset)
    surv_times <- summary(surv_fit_ref)$time
    surv_probs <- summary(surv_fit_ref)$surv
    surv_prob_ref <- numeric(nrow(data_subset))
    for (i in 1:nrow(data_subset)) {
      surv_prob_ref[i] <- approx(surv_times, surv_probs[, i], xout = eval_time, 
                                 method = "linear", rule = 2)$y
    }
    risk_ref <- 1 - surv_prob_ref
    
    # Observed outcomes
    event_ind <- ifelse(data_subset[[event_vars[1]]] <= eval_time & 
                          data_subset[[event_vars[2]]] == 1, 1, 0)
    
    # Brier scores
    brier_new <- mean((surv_prob_new - (1 - event_ind))^2, na.rm = TRUE)
    brier_ref <- mean((surv_prob_ref - (1 - event_ind))^2, na.rm = TRUE)
    
    # Outlier diagnostics
    squared_errors_new <- (surv_prob_new - (1 - event_ind))^2
    outliers <- which(squared_errors_new > 0.9)
    outlier_info <- if (length(outliers) > 0) {
      data.frame(
        Index = data_subset$INDEX[outliers],
        Surv_Prob = round(surv_prob_new[outliers], 4),
        Risk = round(risk_new[outliers], 4),
        Event = event_ind[outliers],
        Squared_Error = round(squared_errors_new[outliers], 4),
        Lin_Pred = round(lin_pred_new[outliers], 4)
      )
    } else {
      NULL
    }
    
    # Prepare data for IDI.INF
    indata <- data.frame(
      time = data_subset[[event_vars[1]]],
      status = data_subset[[event_vars[2]]]
    )
    covs0 <- as.matrix(data_subset[, clin_vars, drop = FALSE])
    covs1 <- as.matrix(data_subset[, c(clin_vars, biomarker_vars), drop = FALSE])
    
    # Check event count
    event_count <- sum(data_subset[[event_vars[2]]] & data_subset[[event_vars[1]]] <= eval_time)
    if (event_count < 5) {
      warning("Too few events (", event_count, ") for reliable IDI.INF computation")
      idi_result <- NULL
    } else {
      # Compute NRI and IDI using IDI.INF
      idi_result <- tryCatch({
        IDI.INF(
          indata = indata,
          covs0 = covs0,
          covs1 = covs1,
          t0 = eval_time,
          npert = 0
        )
      }, error = function(e) {
        warning("Continuous NRI/IDI failed: ", e$message)
        NULL
      })
    }
    
    return(list(
      brier_new = brier_new,
      brier_ref = brier_ref,
      idi_result = idi_result,
      outliers = outlier_info,
      squared_errors_new = squared_errors_new,
      event_count = event_count
    ))
  }
  
  # Compute on full dataset
  full_result <- compute_nri_idi_single(data, new_model_result, ref_formula, eval_time)
  if (is.null(full_result)) {
    stop("NRI/IDI computation failed on full dataset")
  }
  
  # Bootstrap for internal validation
  boot_results <- list(
    brier_new = numeric(n_bootstrap),
    brier_ref = numeric(n_bootstrap),
    nri = numeric(n_bootstrap),
    idi = numeric(n_bootstrap)
  )
  for (i in 1:n_bootstrap) {
    set.seed(i)
    boot_idx <- sample(1:nrow(data), nrow(data), replace = TRUE)
    boot_data <- data[boot_idx, ]
    boot_result <- compute_nri_idi_single(boot_data, new_model_result, ref_formula, eval_time)
    
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
        warning("IDI.INF result incomplete in bootstrap iteration ", i, 
                "; event count: ", boot_result$event_count)
        boot_results$nri[i] <- NA
        boot_results$idi[i] <- NA
      }
    } else {
      warning("Bootstrap iteration ", i, " failed")
      boot_results$brier_new[i] <- NA
      boot_results$brier_ref[i] <- NA
      boot_results$nri[i] <- NA
      boot_results$idi[i] <- NA
    }
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
  
  # Full dataset results
  cat("\nFull Dataset Results:\n")
  cat("Brier Score (New Model):", full_result$brier_new, "\n")
  cat("Brier Score (Reference Model):", full_result$brier_ref, "\n")
  if (!is.null(full_result$idi_result)) {
    cat("\nContinuous NRI and IDI Results:\n")
    print(full_result$idi_result)
  } else {
    cat("\nContinuous NRI and IDI: Not computed\n")
  }
  
  # Outlier diagnostics
  cat("\nOutlier Diagnostics (Full Dataset):\n")
  cat("Summary of squared errors:\n")
  print(summary(full_result$squared_errors_new))
  if (!is.null(full_result$outliers)) {
    cat("Outliers (squared error > 0.9):\n")
    print(full_result$outliers)
  } else {
    cat("No outliers with squared error > 0.9\n")
  }
  
  # Event rate
  cat("\nEvent Rate:\n")
  cat("Events by t=", eval_time, ": ", sum(data$status & data$time <= eval_time), "\n", sep = "")
  cat("Sample size: ", nrow(data), "\n", sep = "")
  
  # Data consistency check
  cat("\nData Consistency Check:\n")
  cat("Summary of lin_pred:\n")
  print(summary(new_model_result$lin_pred))
  cat("Missing values in time/status:", any(is.na(data$time) | is.na(data$status)), "\n")
  cat("Clinical predictors:", clin_vars, "\n")
  cat("Row names alignment check:\n")
  cat("Unique INDEX values in data:", length(unique(data$INDEX)), "\n")
  cat("Unique names in lin_pred:", length(unique(names(new_model_result$lin_pred))), "\n")
  cat("Matching INDEX to lin_pred names:", sum(unique(data$INDEX) %in% names(new_model_result$lin_pred)), "\n")
  
  # Return results
  return(list(
    full_result = full_result,
    boot_results = boot_results
  ))
}
  
# Example usage
# Step 1: Verify clin_vars3 and prot_npx (already defined externally)
# clin_vars3 should be a vector of 3 clinical predictor names, e.g., c("age", "sex", "comorbidity_score")
# prot_npx should be a vector of 21 biomarker names, e.g., c("proteinA", "proteinB", ..., "proteinU")
print(clin_vars3)
print(prot_npx)

# If clin_vars3 shows clinical1 clinical2 clinical3, restore it with your actual predictor names
# Example: clin_vars3 <- c("age", "sex", "comorbidity_score")  # Uncomment and update if needed
# Then verify again:
# print(clin_vars3)

# Step 2: Compute NRI and IDI with bootstrap internal validation
eval_time <- 10
risk_thresholds <- c(0.1, 0.2)  # Example: <10%, 10-20%, >20% risk
nri_idi_result <- compute_nri_idi_internal(
  new_model_result = cox_glmnet_result,  # Use existing cox_glmnet_result
  clin_vars = clin_vars3,
  data = npxdata_all,
  biomarker_vars = prot_npx,
  eval_time = eval_time,
  event_vars = c("time", "status"),
  risk_thresholds = risk_thresholds,
  recalibrate = list(slope = 2.444079, intercept = 1.756214),  # Current calibration
  n_bootstrap = 1000
)

# Step 3: Verify data
cat("\nData Verification:\n")
cat("Summary of lin_pred:\n")
print(summary(cox_glmnet_result$lin_pred))
cat("Missing values in time/status:", any(is.na(npxdata_all$time) | is.na(npxdata_all$status)), "\n")
cat("Predictors in npxdata_all:\n")
print(colnames(npxdata_all))
cat("Clinical predictors:", clin_vars3, "\n")
cat("Biomarkers:", prot_npx, "\n")