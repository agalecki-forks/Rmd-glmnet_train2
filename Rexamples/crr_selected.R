# Define predefined predictors
predefined_vars <- c("var1", "var2", "var3")  # Replace with your actual predictors

# Assume crr_fit is already fitted with predefined_vars
# Example: crr_fit <- crr(npxdata_all$time, npxdata_all$status, npxdata_all[, predefined_vars])

# Create newdata with one row
newdata <- data.frame(matrix(NA, nrow = 1, ncol = length(predefined_vars)))
colnames(newdata) <- predefined_vars

# Calculate mean for numeric variables, mode for factors
if (length(predefined_vars) > 0) {
  for (var in predefined_vars) {
    if (var %in% names(npxdata_all)) {
      if (is.numeric(npxdata_all[[var]])) {
        # Set to mean for numeric variables
        newdata[1, var] <- mean(npxdata_all[[var]], na.rm = TRUE)
      } else if (is.factor(npxdata_all[[var]])) {
        # Set to most frequent level (mode) for factors
        mode_val <- names(sort(table(npxdata_all[[var]]), decreasing = TRUE))[1]
        newdata[1, var] <- factor(mode_val, levels = levels(npxdata_all[[var]]))
      } else {
        warning(paste("Variable", var, "is neither numeric nor factor. Setting to NA."))
        newdata[1, var] <- NA
      }
    } else {
      warning(paste("Variable", var, "not found in npxdata_all. Setting to NA."))
      newdata[1, var] <- NA
    }
  }
}

# Prepare covariates and predict CIF
cov1 <- as.matrix(newdata)
pred <- predict(crr_fit, cov1 = cov1, times = seq(0, max(npxdata_all$time), by = 0.1))

# Plot CIF
plot(pred, xlab = "Time", ylab = "Cumulative Incidence", main = "CIF for CRR Model (Mean/Mode Values)", col = 1)