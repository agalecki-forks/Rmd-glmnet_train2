#' Generate a calibration table with outlier diagnostics
#' @param calib_result Result list from calibration_analysis
#' @param cox_glmnet_result List from cox_glmnet_train
#' @param data Data frame with time and status
#' @param eval_time Time point for calibration
#' @param event_vars Event variable names
#' @return A formatted table and outlier diagnostics
calibration_table_with_diagnostics <- function(calib_result, cox_glmnet_result, data, eval_time, 
                                              event_vars = c("time", "status")) {
  require(survival)
  
  # Extract calibration data
  calib_data <- calib_result$calibration_data
  if (is.null(calib_data) || nrow(calib_data) == 0) {
    stop("No valid calibration data available for table")
  }
  
  # Extract metrics
  calib_slope <- calib_result$calib_slope
  calib_intercept <- calib_result$calib_intercept
  se_slope <- calib_result$se_slope
  se_intercept <- calib_result$se_intercept
  ci_slope <- calib_result$ci_slope
  ci_intercept <- calib_result$ci_intercept
  
  # Calculate individual contributions to Brier score
  bas_cumhazdf <- cox_glmnet_result$bas_cumhazdf
  lin_pred <- cox_glmnet_result$lin_pred
  time <- data[[event_vars[1]]]
  status <- data[[event_vars[2]]]
  
  cumhaz_at_t <- approx(bas_cumhazdf$time, bas_cumhazdf$hazard, xout = eval_time, method = "linear")$y
  surv_prob <- exp(-cumhaz_at_t * exp(lin_pred))
  event_ind <- ifelse(time <= eval_time & status == 1, 1, 0)
  obs_outcome <- ifelse(time > eval_time, 1, event_ind)
  squared_errors <- (surv_prob - obs_outcome)^2
  brier_score <- mean(squared_errors, na.rm = TRUE)
  
  # Calculate bin-level average squared errors
  bins <- cut(surv_prob, breaks = quantile(surv_prob, probs = seq(0, 1, length.out = 11), na.rm = TRUE), 
              include.lowest = TRUE, labels = FALSE)
  bin_errors <- tapply(squared_errors, bins, mean, na.rm = TRUE)
  bin_errors <- bin_errors[match(calib_data$bin, seq_len(10))]
  bin_errors[is.na(bin_errors)] <- 0  # Handle missing bins
  
  # Create calibration table
  calib_table <- data.frame(
    Bin = calib_data$bin,
    Mean_Pred = round(calib_data$mean_pred, 4),
    Obs_Surv = round(calib_data$obs_surv, 4),
    N = calib_data$n,
    Avg_Squared_Error = round(bin_errors, 4)
  )
  
  # Print calibration table
  cat("\nCalibration Table:\n")
  print(calib_table)
  
  # Summary statistics
  cat("\nSummary Statistics:\n")
  cat(sprintf("Calibration Slope: %.3f (SE: %.3f, 95%% CI: %.3f, %.3f)\n", 
              calib_slope, se_slope, ci_slope[1], ci_slope[2]))
  cat(sprintf("Calibration Intercept: %.3f (SE: %.3f, 95%% CI: %.3f, %.3f)\n", 
              calib_intercept, se_intercept, ci_intercept[1], ci_intercept[2]))
  cat(sprintf("Brier Score at t=%d: %.4f\n", eval_time, brier_score))
  
  # Outlier diagnostics
  cat("\nOutlier Diagnostics:\n")
  cat("Summary of squared errors:\n")
  print(summary(squared_errors))
  cat("Top 5 largest squared errors:\n")
  print(head(sort(squared_errors, decreasing = TRUE), 5))
  
  outliers <- which(squared_errors > 0.9)
  if (length(outliers) > 0) {
    cat("Outliers (squared error > 0.9):\n")
    print(data.frame(
      Index = outliers,
      Surv_Prob = round(surv_prob[outliers], 4),
      Obs_Outcome = obs_outcome[outliers],
      Squared_Error = round(squared_errors[outliers], 4),
      Lin_Pred = round(lin_pred[outliers], 4)
    ))
  } else {
    cat("No outliers with squared error > 0.9\n")
  }
  
  # Check lin_pred outliers
  cat("\nLinear Predictor Diagnostics:\n")
  cat("Summary of lin_pred:\n")
  print(summary(lin_pred))
  iqr <- IQR(lin_pred)
  q <- quantile(lin_pred, probs = c(0.25, 0.75))
  extreme_lin_pred <- which(lin_pred < q[1] - 1.5 * iqr | lin_pred > q[2] + 1.5 * iqr)
  if (length(extreme_lin_pred) > 0) {
    cat("Extreme lin_pred values (outside 1.5*IQR):\n")
    print(data.frame(
      Index = extreme_lin_pred,
      Lin_Pred = round(lin_pred[extreme_lin_pred], 4)
    ))
  } else {
    cat("No extreme lin_pred values\n")
  }
  
  # Event rate
  cat("\nEvent Rate:\n")
  cat("Events by t=", eval_time, ": ", sum(data$status & data$time <= eval_time), "\n", sep = "")
  cat("Sample size: ", nrow(data), "\n", sep = "")
  
  # Return table and diagnostics
  return(list(
    calibration_table = calib_table,
    brier_score = brier_score,
    squared_errors = squared_errors,
    bin_errors = bin_errors,
    outliers = outliers
  ))
}

# Example usage with your calibration_analysis output
eval_time <- 10
table_result <- calibration_table_with_diagnostics(
  calib_result = tt, 
  cox_glmnet_result = cox_glmnet_result, 
  data = npxdata_all, 
  eval_time = eval_time, 
  event_vars = c("time", "status")
)

# Print the table again for reference
cat("\nCalibration Table:\n")
print(table_result$calibration_table)
