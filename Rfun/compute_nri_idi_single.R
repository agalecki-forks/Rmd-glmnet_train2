#' Compute NRI and IDI for survival model comparison on a single dataset
#' @param new_model_result List with lin_pred and bas_cumhazdf
#' @param clin_vars Vector of clinical predictor names
#' @param biomarker_vars Vector of biomarker predictor names
#' @param data Data frame with time, status, predictors, and INDEX column
#' @param eval_time Time point for evaluation
#' @param event_vars Event variable names
#' @param recalibrate Apply calibration slope and intercept to new model predictions
#' @return Brier scores, NRI, IDI, diagnostics
compute_nri_idi_single <- function(new_model_result, clin_vars, biomarker_vars, data, eval_time, 
                                   event_vars = c("time", "status"), 
                                   recalibrate = list(slope = 1, intercept = 0)) {
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
  
  # Standardize covariates
  data_scaled <- data
  for (var in all_vars) {
    data_scaled[[var]] <- scale(data_scaled[[var]], center = TRUE, scale = TRUE)
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
  
  # Align lin_pred using INDEX column
  lin_pred_new <- new_model_result$lin_pred[match(as.character(data_scaled$INDEX), 
                                                  as.character(names(new_model_result$lin_pred)))]
  if (any(is.na(lin_pred_new))) {
    stop("Missing lin_pred values for ", sum(is.na(lin_pred_new)), " observations.")
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
  ref_formula <- as.formula(paste("Surv(", event_vars[1], ",", event_vars[2], ") ~", 
                                  paste(clin_vars, collapse = " + ")))
  ref_model <- tryCatch({
    environment(ref_formula) <- list2env(list(data_subset = data_scaled), 
                                         parent = .GlobalEnv)
    coxph(ref_formula, data = data_scaled, x = TRUE, 
          control = coxph.control(iter.max = 100))
  }, error = function(e) {
    warning("Reference model fitting failed: ", e$message)
    return(NULL)
  })
  if (is.null(ref_model)) {
    return(NULL)
  }
  
  # Predict survival probabilities for reference model
  surv_fit_ref <- survfit(ref_model, newdata = data_scaled)
  surv_times <- summary(surv_fit_ref)$time
  surv_probs <- summary(surv_fit_ref)$surv
  surv_prob_ref <- numeric(nrow(data_scaled))
  for (i in 1:nrow(data_scaled)) {
    surv_prob_ref[i] <- approx(surv_times, surv_probs[, i], xout = eval_time, 
                               method = "linear", rule = 2)$y
  }
  risk_ref <- 1 - surv_prob_ref
  
  # Observed outcomes
  event_ind <- ifelse(data_scaled[[event_vars[1]]] <= eval_time & 
                        data_scaled[[event_vars[2]]] == 1, 1, 0)
  
  # Brier scores
  brier_new <- mean((surv_prob_new - (1 - event_ind))^2, na.rm = TRUE)
  brier_ref <- mean((surv_prob_ref - (1 - event_ind))^2, na.rm = TRUE)
  
  # Outlier diagnostics
  squared_errors_new <- (surv_prob_new - (1 - event_ind))^2
  outliers <- which(squared_errors_new > 0.9)
  outlier_info <- if (length(outliers) > 0) {
    data.frame(
      Index = data_scaled$INDEX[outliers],
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
    time = data_scaled[[event_vars[1]]],
    status = data_scaled[[event_vars[2]]]
  )
  covs0 <- as.matrix(data_scaled[, clin_vars, drop = FALSE])
  covs1 <- as.matrix(data_scaled[, c(clin_vars, biomarker_vars), drop = FALSE])
  
  # Check event count
  event_count <- sum(data_scaled[[event_vars[2]]] & data_scaled[[event_vars[1]]] <= eval_time)
  if (event_count < 10) {
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
  
  # Print results
  cat("\nSingle Dataset Results:\n")
  cat("Brier Score (New Model):", brier_new, "\n")
  cat("Brier Score (Reference Model):", brier_ref, "\n")
  if (!is.null(idi_result)) {
    cat("\nContinuous NRI and IDI Results:\n")
    #cat("IDI:", idi_result$point$IDI[1,1], "\n")
    #cat("NRI:", idi_result$point$TXI[1], "\n")
    # Removed: print(idi_result)
    print(str(idi_result))
  } else {
    cat("\nContinuous NRI and IDI: Not computed\n")
  }
  
  cat("\nOutlier Diagnostics:\n")
  cat("Summary of squared errors:\n")
  print(summary(squared_errors_new))
  if (!is.null(outlier_info)) {
    cat("Outliers (squared error > 0.9):\n")
    print(outlier_info)
  } else {
    cat("No outliers with squared error > 0.9\n")
  }
  
  cat("\nEvent Rate:\n")
  cat("Events by t=", eval_time, ": ", event_count, "\n", sep = "")
  cat("Sample size: ", nrow(data), "\n", sep = "")
  
  cat("\nData Consistency Check:\n")
  cat("Summary of lin_pred:\n")
  print(summary(new_model_result$lin_pred))
  cat("Missing values in time/status:", any(is.na(data$time) | is.na(data$status)), "\n")
  cat("Clinical predictors:", clin_vars, "\n")
  cat("Row names alignment check:\n")
  cat("Unique INDEX values in data:", length(unique(data$INDEX)), "\n")
  cat("Unique names in lin_pred:", length(unique(names(new_model_result$lin_pred))), "\n")
  cat("Matching INDEX to lin_pred names:", sum(unique(data$INDEX) %in% names(new_model_result$lin_pred)), "\n")
  
  return(list(
    brier_new = brier_new,
    brier_ref = brier_ref,
    idi_result = idi_result,
    outliers = outlier_info,
    squared_errors_new = squared_errors_new,
    event_count = event_count
  ))
}

# Example:
#eval_time <- 10
#single_result <- compute_nri_idi_single(
#  new_model_result = cox_glmnet_result,
#  clin_vars = clin_vars3,
#  biomarker_vars = prot_npx, # top_biomarkers,
#  data = npxdata_all,
#  eval_time = eval_time,
#  event_vars = c("time", "status"),
#  recalibrate = list(slope = 2.444079, intercept = 1.756214)
# )
#cat("IDI:", idi_result$point$IDI[1,1], "\n")
#cat("NRI:", idi_result$point$TXI[1], "\n")
# Removed: print(idi_result)
