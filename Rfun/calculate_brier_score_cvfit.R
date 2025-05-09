#' Calculate Brier score for a penalized Cox model using cox_glmnet_train output
#' @param cox_glmnet_result List returned by cox_glmnet_train (contains cvfit_scaled, bas_cumhazdf, lin_pred, scaling_params)
#' @param data Data frame containing time and status variables
#' @param eval_time Numeric, time point at which to evaluate the Brier score (e.g., 10)
#' @param event_vars Character vector of event variables (default: c("time", "status"))
#' @return Numeric, Brier score at the specified evaluation time
calculate_brier_score <- function(cox_glmnet_result, data, eval_time, event_vars = c("time", "status")) {
  # Load required package
  require(survival)
  
  # Extract components from cox_glmnet_result
  bas_cumhazdf <- cox_glmnet_result$bas_cumhazdf
  lin_pred <- cox_glmnet_result$lin_pred
  
  # Validate inputs
  if (!all(event_vars %in% colnames(data))) {
    stop("Event variables not found in data")
  }
  
  # Extract time and status
  time <- data[[event_vars[1]]]
  status <- data[[event_vars[2]]]
  
  # Diagnostic: Check linear predictors
  if (any(is.na(lin_pred)) || any(is.nan(lin_pred))) {
    stop("Linear predictors contain NA or NaN values.")
  }
  if (max(abs(lin_pred)) > 10) {
    warning("Linear predictors have extreme values (range: ", min(lin_pred), " to ", max(lin_pred), 
            "). This may lead to unstable survival probabilities.")
  }
  
  # Get baseline cumulative hazard at eval_time
  haz_times <- bas_cumhazdf$time
  haz_values <- bas_cumhazdf$hazard
  if (eval_time > max(haz_times)) {
    warning("eval_time exceeds maximum time in bas_cumhazdf; using maximum available time")
    eval_time <- max(haz_times)
  }
  cumhaz_at_t <- approx(haz_times, haz_values, xout = eval_time, method = "linear")$y
  
  # Calculate survival probability: S(t) = exp(-H(t) * exp(lin_pred))
  surv_prob <- exp(-cumhaz_at_t * exp(lin_pred))
  
  # Diagnostic: Check survival probabilities
  if (all(is.na(surv_prob)) || any(is.nan(surv_prob))) {
    stop("Survival probabilities contain NA or NaN values. Check baseline hazard or linear predictors.")
  }
  if (length(unique(surv_prob)) == 1) {
    warning("All survival probabilities are identical (", unique(surv_prob), 
            "). Brier score may reflect poor model discrimination.")
  }
  
  # Create indicator for event or censoring
  event_ind <- ifelse(time <= eval_time & status == 1, 1, 0)
  obs_outcome <- ifelse(time > eval_time, 1, event_ind)
  
  # Calculate Brier score
  brier_score <- mean((surv_prob - obs_outcome)^2, na.rm = TRUE)
  
  return(brier_score)
}

# -- Example usage
#eval_time <- 10  # Evaluate Brier score at t=10
#brier_score <- calculate_brier_score(cox_glmnet_result = cox_glmnet_result, 
#                                     data = npxdata_all, 
#                                     eval_time = eval_time, 
#                                     event_vars = c("time", "status"))
#cat("Brier score at time", eval_time, ":", brier_score, "\n")

# Diagnostic output
#lin_pred <- cox_glmnet_result$lin_pred
#surv_prob <- exp(-approx(cox_glmnet_result$bas_cumhazdf$time, 
#                         cox_glmnet_result$bas_cumhazdf$hazard, 
#                         xout = eval_time, method = "linear")$y * exp(lin_pred))
#cat("Range of lin_pred:", range(lin_pred, na.rm = TRUE), "\n")
#cat("Number of unique survival probabilities:", length(unique(surv_prob)), "\n")
# summary(surv_prob)