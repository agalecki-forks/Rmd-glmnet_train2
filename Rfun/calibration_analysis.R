#' Calibration analysis using cox_glmnet_train output
#' @param cox_glmnet_result List from cox_glmnet_train
#' @param data Data frame with time and status
#' @param eval_time Time point for calibration
#' @param event_vars Event variable names
#' @param n_bins Number of bins
#' @param clip_value Clipping value for boundary adjustment
#' @param exclude_extremes Exclude bins with obs_surv = 0 or 1
#' @return List with calibration data, slope, intercept, standard errors, and CIs
calibration_analysis <- function(cox_glmnet_result, data, eval_time, event_vars = c("time", "status"), 
                                n_bins = 10, clip_value = 1e-3, exclude_extremes = TRUE) {
  require(survival)
  
  # Extract components
  bas_cumhazdf <- cox_glmnet_result$bas_cumhazdf
  lin_pred <- cox_glmnet_result$lin_pred
  
  # Validate inputs
  if (!all(event_vars %in% colnames(data))) {
    stop("Event variables not found in data")
  }
  
  time <- data[[event_vars[1]]]
  status <- data[[event_vars[2]]]
  
  # Diagnostic: Check linear predictors
  if (any(is.na(lin_pred)) || any(is.nan(lin_pred))) {
    stop("Linear predictors contain NA or NaN values.")
  }
  if (max(abs(lin_pred)) > 10) {
    warning("Extreme lin_pred range: ", min(lin_pred), " to ", max(lin_pred))
  }
  
  # Interpolate baseline hazard
  haz_times <- bas_cumhazdf$time
  haz_values <- bas_cumhazdf$hazard
  if (eval_time > max(haz_times)) {
    warning("eval_time exceeds maximum time in bas_cumhazdf; using maximum available time")
    eval_time <- max(haz_times)
  }
  cumhaz_at_t <- approx(haz_times, haz_values, xout = eval_time, method = "linear")$y
  
  # Calculate survival probabilities
  pred_surv_prob <- as.numeric(exp(-cumhaz_at_t * exp(lin_pred)))
  
  # Diagnostic: Check survival probabilities
  if (!is.numeric(pred_surv_prob)) {
    stop("pred_surv_prob is not numeric")
  }
  if (any(is.na(pred_surv_prob)) || any(is.nan(pred_surv_prob))) {
    stop("Survival probabilities contain NA or NaN values.")
  }
  if (length(unique(pred_surv_prob)) == 1) {
    warning("All survival probabilities identical (", unique(pred_surv_prob), "). Calibration cannot proceed.")
    return(list(calibration_data = NULL, calib_slope = NA, calib_intercept = NA, 
                se_slope = NA, se_intercept = NA, ci_slope = c(NA, NA), ci_intercept = c(NA, NA)))
  }
  
  # Create calibration data frame
  calib_data <- data.frame(
    pred_surv_prob = as.numeric(pred_surv_prob),
    time = time,
    status = status,
    stringsAsFactors = FALSE
  )
  
  # Bin predicted probabilities
  breaks <- tryCatch({
    quantile(pred_surv_prob, probs = seq(0, 1, length.out = n_bins + 1), na.rm = TRUE)
  }, error = function(e) {
    warning("Quantile binning failed; using equal-width bins: ", e$message)
    seq(min(pred_surv_prob, na.rm = TRUE), max(pred_surv_prob, na.rm = TRUE), length.out = n_bins + 1)
  })
  
  unique_breaks <- unique(breaks)
  if (length(unique_breaks) < 2) {
    warning("Fewer than 2 unique breaks; calibration cannot proceed.")
    return(list(calibration_data = NULL, calib_slope = NA, calib_intercept = NA, 
                se_slope = NA, se_intercept = NA, ci_slope = c(NA, NA), ci_intercept = c(NA, NA)))
  }
  if (length(unique_breaks) < n_bins + 1) {
    warning(sprintf("Only %d unique breaks; reducing to %d bins", length(unique_breaks) - 1, length(unique_breaks) - 1))
    n_bins <- length(unique_breaks) - 1
    breaks <- unique_breaks
  }
  
  # Bin predicted probabilities
  calib_data$bin <- tryCatch({
    cut(pred_surv_prob, breaks = breaks, include.lowest = TRUE, labels = FALSE)
  }, error = function(e) {
    warning("Binning failed; calibration cannot proceed: ", e$message)
    return(list(calibration_data = NULL, calib_slope = NA, calib_intercept = NA, 
                se_slope = NA, se_intercept = NA, ci_slope = c(NA, NA), ci_intercept = c(NA, NA)))
  })
  
  if (is.null(calib_data$bin)) {
    return(list(calibration_data = NULL, calib_slope = NA, calib_intercept = NA, 
                se_slope = NA, se_intercept = NA, ci_slope = c(NA, NA), ci_intercept = c(NA, NA)))
  }
  
  # Diagnostic: Check binning
  cat("Bin counts:\n")
  print(table(calib_data$bin, useNA = "always"))
  
  # Calculate observed survival probabilities per bin
  calib_summary <- do.call(rbind, lapply(1:n_bins, function(b) {
    bin_data <- calib_data[calib_data$bin == b & !is.na(calib_data$bin), , drop = FALSE]
    if (nrow(bin_data) == 0) {
      warning("Bin ", b, " is empty")
      return(data.frame(bin = b, mean_pred = NA, obs_surv = NA, n = 0))
    }
    if (!is.numeric(bin_data$pred_surv_prob)) {
      warning("Non-numeric pred_surv_prob in bin ", b, ": ", class(bin_data$pred_surv_prob))
      return(data.frame(bin = b, mean_pred = NA, obs_surv = NA, n = nrow(bin_data)))
    }
    if (any(is.na(bin_data$pred_surv_prob))) {
      warning("NA values in pred_surv_prob in bin ", b)
    }
    km_fit <- tryCatch({
      survfit(Surv(time, status) ~ 1, data = bin_data)
    }, error = function(e) {
      warning("survfit failed for bin ", b, ": ", e$message)
      return(NULL)
    })
    if (is.null(km_fit)) {
      return(data.frame(bin = b, mean_pred = NA, obs_surv = NA, n = nrow(bin_data)))
    }
    km_summary <- summary(km_fit, times = eval_time)
    obs_surv <- ifelse(length(km_summary$surv) > 0 && !is.na(km_summary$surv), km_summary$surv, NA)
    data.frame(
      bin = b,
      mean_pred = mean(bin_data$pred_surv_prob, na.rm = TRUE),
      obs_surv = obs_surv,
      n = nrow(bin_data)
    )
  }))
  
  # Remove NA rows
  calib_summary <- calib_summary[!is.na(calib_summary$obs_surv) & !is.na(calib_summary$mean_pred), ]
  
  # Calibration plot
  if (nrow(calib_summary) > 0) {
    plot(calib_summary$mean_pred, calib_summary$obs_surv, 
         xlab = "Mean Predicted Survival Probability", 
         ylab = "Observed Survival Probability (Kaplan-Meier)",
         main = paste("Calibration Plot at Time", eval_time),
         xlim = c(0, 1), ylim = c(0, 1), pch = 19, col = "blue")
    abline(0, 1, lty = 2, col = "red")
    text(calib_summary$mean_pred, calib_summary$obs_surv, labels = calib_summary$n, pos = 3, cex = 0.8)
  } else {
    warning("No valid bins for calibration plot")
    return(list(calibration_data = calib_summary, calib_slope = NA, calib_intercept = NA, 
                se_slope = NA, se_intercept = NA, ci_slope = c(NA, NA), ci_intercept = c(NA, NA)))
  }
  
  # Calibration slope and intercept
  calib_slope <- NA
  calib_intercept <- NA
  se_slope <- NA
  se_intercept <- NA
  ci_slope <- c(NA, NA)
  ci_intercept <- c(NA, NA)
  
  if (nrow(calib_summary) > 1) {
    # Adjust mean_pred and obs_surv to avoid boundary values
    calib_summary$mean_pred_adj <- pmin(pmax(calib_summary$mean_pred, clip_value), 1 - clip_value)
    calib_summary$obs_surv_adj <- pmin(pmax(calib_summary$obs_surv, clip_value), 1 - clip_value)
    
    calib_summary$logit_pred <- log(calib_summary$mean_pred_adj / (1 - calib_summary$mean_pred_adj))
    calib_summary$logit_obs <- log(calib_summary$obs_surv_adj / (1 - calib_summary$obs_surv_adj))
    
    # Diagnostic: Check logit values
    cat("Logit values:\n")
    print(calib_summary[, c("bin", "mean_pred", "obs_surv", "logit_pred", "logit_obs")])
    
    # Optionally exclude bins with obs_surv = 0 or 1
    calib_summary_reg <- calib_summary
    if (exclude_extremes) {
      calib_summary_reg <- calib_summary[calib_summary$obs_surv > 0 & calib_summary$obs_surv < 1, ]
      if (nrow(calib_summary_reg) <= 1) {
        warning("Insufficient valid bins after excluding extremes")
        return(list(calibration_data = calib_summary, calib_slope = NA, calib_intercept = NA, 
                    se_slope = NA, se_intercept = NA, ci_slope = c(NA, NA), ci_intercept = c(NA, NA)))
      }
    }
    
    if (all(is.finite(calib_summary_reg$logit_pred)) && all(is.finite(calib_summary_reg$logit_obs))) {
      calib_model <- lm(logit_obs ~ logit_pred, data = calib_summary_reg, weights = n)
      calib_slope <- coef(calib_model)[2]
      calib_intercept <- coef(calib_model)[1]
      # Extract standard errors
      model_summary <- summary(calib_model)
      se_intercept <- model_summary$coefficients[1, "Std. Error"]
      se_slope <- model_summary$coefficients[2, "Std. Error"]
      # Compute 95% confidence intervals
      ci_slope <- calib_slope + c(-1.96, 1.96) * se_slope
      ci_intercept <- calib_intercept + c(-1.96, 1.96) * se_intercept
    } else {
      warning("Non-finite logit values after adjustment; cannot compute calibration slope/intercept")
    }
  } else {
    warning("Insufficient valid bins for calibration slope/intercept")
  }
  
  return(list(calibration_data = calib_summary, 
              calib_slope = calib_slope, 
              calib_intercept = calib_intercept, 
              se_slope = se_slope, 
              se_intercept = se_intercept,
              ci_slope = ci_slope,
              ci_intercept = ci_intercept))
}

#--- Example usage
#eval_time <- 10
#tt <- calibration_analysis(cox_glmnet_result = cox_glmnet_result, 
#                           data = npxdata_all, 
#                           eval_time = eval_time, 
#                           event_vars = c("time", "status"), 
#                           n_bins = 10, 
#                           clip_value = 1e-3, 
#                           exclude_extremes = TRUE)
#cat("Calibration Slope:", tt$calib_slope, "\n")
#cat("Calibration Intercept:", tt$calib_intercept, "\n")
#cat("Standard Error (Slope):", tt$se_slope, "\n")
#cat("Standard Error (Intercept):", tt$se_intercept, "\n")
#cat("Slope 95% CI:", tt$ci_slope[1], "to", tt$ci_slope[2], "\n")
#cat("Intercept 95% CI:", tt$ci_intercept[1], "to", tt$ci_intercept[2], "\n")
#if (!is.null(tt$calibration_data)) {
#  print(tt$calibration_data)
#}
# plot(calib_data$mean_pred, calib_data$obs_surv ...