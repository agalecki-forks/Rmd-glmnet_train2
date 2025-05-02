# Function to check survival data consistency with flexible variable names
check_survival_consistency <- function(data, 
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
  
  # Check for invalid values
  issues <- list()
  
  # 1. Check if status and death are 0/1
  if (any(!temp_data$status %in% c(0, 1), na.rm = TRUE)) {
    issues <- c(issues, paste(status_var, "contains values other than 0 or 1"))
  }
  if (any(!temp_data$death %in% c(0, 1), na.rm = TRUE)) {
    issues <- c(issues, paste(death_var, "contains values other than 0 or 1"))
  }
  
  # 2. Check if time and time_death are non-negative
  if (any(temp_data$time < 0, na.rm = TRUE)) {
    issues <- c(issues, paste(time_var, "contains negative values"))
  }
  if (any(temp_data$time_death < 0, na.rm = TRUE)) {
    issues <- c(issues, paste(time_death_var, "contains negative values"))
  }
  
  # 3. Check consistency: time_death >= time when death == 1
  invalid_time <- temp_data$death == 1 & temp_data$time_death < temp_data$time
  if (any(invalid_time, na.rm = TRUE)) {
    invalid_indices <- temp_data$INDEX[invalid_time]
    issues <- c(issues, paste("Death before time for", index_var, ":", 
                              paste(invalid_indices, collapse = ", ")))
  }
  
  # 4. Check if status aligns with death
  status_mismatch <- (temp_data$status == 1 & temp_data$death == 0) | 
                     (temp_data$status == 0 & temp_data$death == 1)
  if (any(status_mismatch, na.rm = TRUE)) {
    mismatch_indices <- temp_data$INDEX[status_mismatch]
    issues <- c(issues, paste(status_var, "and", death_var, "mismatch for", index_var, ":", 
                              paste(mismatch_indices, collapse = ", ")))
  }
  
  # Return results
  if (length(issues) == 0) {
    message("All checks passed: Data is consistent.")
  } else {
    warning("Issues found:\n", paste("-", issues, collapse = "\n"))
  }
}
