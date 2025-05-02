summarize_data <- function(npx_data) {
  # Validate input
  if (!is.data.frame(npx_data)) stop("Input 'npx_data' must be a data frame")
  
  # Get numeric variables
  numeric_vars <- names(npx_data)[sapply(npx_data, is.numeric)]
  if (length(numeric_vars) == 0) stop("No numeric variables found in 'npx_data'")
  
  # Check for problematic variable names (e.g., "25%")
  if (any(grepl("^25%", numeric_vars))) {
    warning("Variable names starting with '25%' detected: ", 
            paste(numeric_vars[grepl("^25%", numeric_vars)], collapse = ", "))
  }
  
  # Compute summary statistics
  summary_stats <- data.frame(
    do.call(rbind, lapply(numeric_vars, function(var) {
      summarize_variable(npx_data[[var]], var)
    })),
    row.names = NULL
  )
  
  # Round only appropriate columns (mean, SD, min, pct25, pct50, pct75, max)
  cols_to_round <- c("mean", "SD", "min", "pct25", "pct50", "pct75", "max")
  summary_stats[, cols_to_round] <- lapply(summary_stats[, cols_to_round], round, digits = 2)
  
  return(summary_stats)
}
