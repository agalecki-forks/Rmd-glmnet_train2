# Function to check survival data consistency with detailed output
# Function to check survival data consistency for two events (e.g., ESKD and death)
check_survival_consistency2 <- function(data, 
                                       index_var = "INDEX", 
                                       time_var = "time", 
                                       status_var = "status", 
                                       time_death_var = "time_death", 
                                       death_var = "death") {
  # Map provided variable names to expected ones
  required_vars <- c(index_var, time_var, status_var, time_death_var, death_var)
  missing_vars <- setdiff(required_vars, names(data))
  if (length(missing_vars) > 0) {
    stop("Missing variables: ", paste(missing_vars, collapse = ", "))
  }
  
  # Create a temporary dataset with standardized names
  temp_data <- data
  names(temp_data)[names(temp_data) == index_var] <- "INDEX"
  names(temp_data)[names(temp_data) == time_var] <- "time"
  names(temp_data)[names(temp_data) == status_var] <- "status"
  names(temp_data)[names(temp_data) == time_death_var] <- "time_death"
  names(temp_data)[names(temp_data) == death_var] <- "death"
  
  # Initialize results list
  results <- list(
    summary = list(),
    issues = list(),
    problematic_rows = data.frame()
  )
  
  # Summary statistics
  results$summary$n_rows <- nrow(temp_data)
  results$summary$missing_time <- sum(is.na(temp_data$time))
  results$summary$missing_time_death <- sum(is.na(temp_data$time_death))
  results$summary$missing_status <- sum(is.na(temp_data$status))
  results$summary$missing_death <- sum(is.na(temp_data$death))
  results$summary$first_event_count <- sum(temp_data$status == 1, na.rm = TRUE)
  results$summary$death_event_count <- sum(temp_data$death == 1, na.rm = TRUE)
  results$summary$censored <- sum(temp_data$status == 0, na.rm = TRUE)
  
  # Check for invalid values
  # 1. Check if status and death are 0/1
  invalid_status <- !temp_data$status %in% c(0, 1) & !is.na(temp_data$status)
  if (any(invalid_status)) {
    n_invalid <- sum(invalid_status)
    results$issues <- c(results$issues, 
                        paste(status_var, "has", n_invalid, "non-0/1 values"))
    results$problematic_rows <- rbind(
      results$problematic_rows,
      data.frame(
        Issue = paste(status_var, "non-0/1"),
        INDEX = temp_data$INDEX[invalid_status],
        Value = temp_data$status[invalid_status]
      )
    )
  }
  
  invalid_death <- !temp_data$death %in% c(0, 1) & !is.na(temp_data$death)
  if (any(invalid_death)) {
    n_invalid <- sum(invalid_death)
    results$issues <- c(results$issues, 
                        paste(death_var, "has", n_invalid, "non-0/1 values"))
    results$problematic_rows <- rbind(
      results$problematic_rows,
      data.frame(
        Issue = paste(death_var, "non-0/1"),
        INDEX = temp_data$INDEX[invalid_death],
        Value = temp_data$death[invalid_death]
      )
    )
  }
  
  # 2. Check if time and time_death are non-negative
  negative_time <- temp_data$time < 0 & !is.na(temp_data$time)
  if (any(negative_time)) {
    n_negative <- sum(negative_time)
    results$issues <- c(results$issues, 
                        paste(time_var, "has", n_negative, "negative values"))
    results$problematic_rows <- rbind(
      results$problematic_rows,
      data.frame(
        Issue = paste(time_var, "negative"),
        INDEX = temp_data$INDEX[negative_time],
        Value = temp_data$time[negative_time]
      )
    )
  }
  
  negative_time_death <- temp_data$time_death < 0 & !is.na(temp_data$time_death)
  if (any(negative_time_death)) {
    n_negative <- sum(negative_time_death)
    results$issues <- c(results$issues, 
                        paste(time_death_var, "has", n_negative, "negative values"))
    results$problematic_rows <- rbind(
      results$problematic_rows,
      data.frame(
        Issue = paste(time_death_var, "negative"),
        INDEX = temp_data$INDEX[negative_time_death],
        Value = temp_data$time_death[negative_time_death]
      )
    )
  }
  
  # 3. Check consistency: time_death >= time when death == 1
  invalid_time <- temp_data$death == 1 & temp_data$time_death < temp_data$time & 
                 !is.na(temp_data$time) & !is.na(temp_data$time_death)
  if (any(invalid_time)) {
    n_invalid <- sum(invalid_time)
    results$issues <- c(results$issues, 
                        paste(n_invalid, "cases where", time_death_var, "<", 
                              time_var, "when", death_var, "= 1"))
    results$problematic_rows <- rbind(
      results$problematic_rows,
      data.frame(
        Issue = paste(time_death_var, "<", time_var),
        INDEX = temp_data$INDEX[invalid_time],
        Value = paste(temp_data$time_death[invalid_time], "<", 
                      temp_data$time[invalid_time])
      )
    )
  }
  
  # Print detailed output
  cat("=== Survival Data Consistency Check (Two Events) ===\n")
  cat("Dataset:", deparse(substitute(data)), "\n")
  cat("Variable Names:\n")
  cat("  Index:", index_var, "\n")
  cat("  Time (First Event or Censoring):", time_var, "\n")
  cat("  Status (First Event, e.g., ESKD):", status_var, "\n")
  cat("  Time of Death:", time_death_var, "\n")
  cat("  Death Event:", death_var, "\n\n")
  
  cat("Summary Statistics:\n")
  cat("  Total Rows:", results$summary$n_rows, "\n")
  cat("  Missing", time_var, ":", results$summary$missing_time, "\n")
  cat("  Missing", time_death_var, ":", results$summary$missing_time_death, "\n")
  cat("  Missing", status_var, ":", results$summary$missing_status, "\n")
  cat("  Missing", death_var, ":", results$summary$missing_death, "\n")
  cat("  First Event (e.g., ESKD) (", status_var, "=1):", 
      results$summary$first_event_count, "\n")
  cat("  Death Events (", death_var, "=1):", 
      results$summary$death_event_count, "\n")
  cat("  Censored (", status_var, "=0):", results$summary$censored, "\n\n")
  
  if (length(results$issues) == 0) {
    cat("Result: All checks passed! Data is consistent.\n")
    cat("Note: No check for alignment between", status_var, "and", death_var, 
        "as they represent different events (e.g., ESKD and death).\n")
  } else {
    cat("Result: Issues detected!\n")
    cat("Issues Found:\n")
    for (i in seq_along(results$issues)) {
      cat("  ", i, ". ", results$issues[[i]], "\n")
    }
    cat("\nProblematic Rows (First 5 per issue, if any):\n")
    if (nrow(results$problematic_rows) > 0) {
      # Group by issue and show up to 5 rows per issue
      issues <- unique(results$problematic_rows$Issue)
      for (issue in issues) {
        issue_rows <- results$problematic_rows[results$problematic_rows$Issue == issue, ]
        cat("  Issue:", issue, "\n")
        print(head(issue_rows[, c("INDEX", "Value")], 5))
        cat("\n")
      }
    } else {
      cat("  None\n")
    }
 #   cat("Note: No check for alignment between", status_var, "and", death_var, 
 #       "as they represent different events (e.g., ESKD and death).\n")
  }
  
  # Return results invisibly for further use if needed
  invisible(results)
}
