calculate_survival_probabilities <- function(lin_pred, baseline_df) {
  # Ensure baseline_df is sorted by time
  baseline_df <- baseline_df[order(baseline_df$time), ]
  
  # Compute cumulative baseline hazard (H_0(t))
  time_diff <- c(baseline_df$time[1], diff(baseline_df$time))
  cumulative_hazard <- cumsum(baseline_df$baseline_hazard * time_diff)
  
  # Compute baseline survival function S_0(t) = exp(-H_0(t))
  baseline_survival <- exp(-cumulative_hazard)
  
  # Initialize matrix for survival probabilities
  survival_probs <- matrix(NA, nrow = length(baseline_df$time), ncol = length(lin_pred))
  rownames(survival_probs) <- baseline_df$time
  colnames(survival_probs) <- paste0("individual_", 1:length(lin_pred))
  
  # Compute survival probabilities for each individual
  for (i in 1:length(lin_pred)) {
    survival_probs[, i] <- baseline_survival ^ exp(lin_pred[i])
  }
  
  # Return as data frame
  return(as.data.frame(survival_probs))
}

