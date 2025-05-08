# Load required packages
library(survival)
library(pROC)
library(survivalROC)
library(glmnet)
library(dplyr)

# Set seed for reproducibility
set.seed(123)

# Assume npxdata_all is already loaded with columns: time, status, clin_vars, prot_npx
# Example structure:
# npxdata_all <- data.frame(
#   time = ..., status = ..., 
#   age = ..., sex = ..., bmi = ...,  # clin_vars
#   prot1 = ..., prot2 = ..., ..., prot21 = ...  # prot_npx
# )
# clin_vars <- c("age", "sex", "bmi")
# prot_npx <- paste0("prot", 1:21)

# Extract selected proteins from cv.glmnet object (cvfit)
# Get coefficients at lambda.min
coef_selected <- coef(cvfit, s = "lambda.min")
selected_proteins <- rownames(coef_selected)[which(coef_selected != 0)]
cat("Selected proteins:", selected_proteins, "\n")

# Ensure selected_proteins is limited to those in prot_npx
selected_proteins <- intersect(selected_proteins, prot_npx)
if (length(selected_proteins) != 8) {
  warning("Number of selected proteins is not 8. Check cvfit or prot_npx.")
}

# Prepare data for modeling
# Ensure no missing data
npxdata_all <- npxdata_all %>% 
  filter(complete.cases(.))  # Remove rows with NA

# Define survival outcome
surv_obj <- Surv(npxdata_all$time, npxdata_all$status)

# --- Model 1: Combined model (3 clinical predictors + 8 selected proteins) ---
# Create formula for combined model
combined_vars <- c(clin_vars, selected_proteins)
combined_formula <- as.formula(
  paste("surv_obj ~", paste(combined_vars, collapse = " + "))
)

# Fit Cox model
cox_combined <- coxph(combined_formula, data = npxdata_all)

# Predict linear predictor (risk score)
pred_combined <- predict(cox_combined, type = "lp")

# --- Models 2-22: Individual protein models (3 clinical predictors + 1 protein) ---
# Initialize list to store models and predictions
individual_models <- list()
pred_individual <- list()

# Loop over each protein in prot_npx
for (protein in prot_npx) {
  # Create formula
  ind_formula <- as.formula(
    paste("surv_obj ~", paste(c(clin_vars, protein), collapse = " + "))
  )
  
  # Fit Cox model
  cox_ind <- coxph(ind_formula, data = npxdata_all)
  
  # Store model
  individual_models[[protein]] <- cox_ind
  
  # Predict linear predictor
  pred_individual[[protein]] <- predict(cox_ind, type = "lp")
}

# --- Compare Predictive Power ---
# Use time-dependent AUC at a specific time point (e.g., median survival time)
time_point <- median(npxdata_all$time, na.rm = TRUE)
cat("Time point for AUC:", time_point, "\n")

# Function to compute time-dependent AUC using survivalROC
compute_auc <- function(predictions, time, status, time_point) {
  roc <- survivalROC(
    Stime = time,
    status = status,
    marker = predictions,
    predict.time = time_point,
    method = "NNE"  # Nearest Neighbor Estimation
  )
  return(roc$AUC)
}

# Compute AUC for combined model
auc_combined <- compute_auc(
  pred_combined, npxdata_all$time, npxdata_all$status, time_point
)

# Compute AUC for each individual model
auc_individual <- sapply(names(pred_individual), function(protein) {
  compute_auc(
    pred_individual[[protein]], npxdata_all$time, npxdata_all$status, time_point
  )
})

# --- Perform DeLong’s Test for AUC Comparisons ---
# Create ROC objects for DeLong’s test (using risk scores and status at time_point)
# Subset to subjects still at risk at time_point
at_risk <- npxdata_all$time >= time_point | npxdata_all$status == 1
status_at_time <- npxdata_all$status[at_risk]

# Combined model ROC
roc_combined <- roc(
  status_at_time,
  pred_combined[at_risk],
  quiet = TRUE
)

# Individual model ROCs and DeLong tests
delong_results <- lapply(names(pred_individual), function(protein) {
  roc_ind <- roc(
    status_at_time,
    pred_individual[[protein]][at_risk],
    quiet = TRUE
  )
  test <- roc.test(roc_combined, roc_ind, method = "delong")
  return(list(
    protein = protein,
    auc_ind = auc_individual[protein],
    p_value = test$p.value,
    z_stat = test$statistic
  ))
})

# Convert DeLong results to data frame
delong_df <- do.call(rbind, lapply(delong_results, as.data.frame))
delong_df$auc_combined <- auc_combined

# --- Harrell’s Concordance (Optional) ---
# Concordance for combined model
conc_combined <- concordance(cox_combined)$concordance

# Concordance for individual models
conc_individual <- sapply(individual_models, function(model) {
  concordance(model)$concordance
})

# --- Summarize Results ---
# Create results table
results_table <- data.frame(
  Protein = c("Combined", names(auc_individual)),
  AUC = c(auc_combined, auc_individual),
  Concordance = c(conc_combined, conc_individual),
  DeLong_p_value = c(NA, delong_df$p_value),
  DeLong_Z = c(NA, delong_df$z_stat)
)

# Round numeric columns
results_table[, 2:5] <- lapply(results_table[, 2:5], round, digits = 3)

# Print results
print("Comparison of Predictive Performance:")
print(results_table)

# --- Visualize Results ---
# Bar plot of AUCs
pdf("auc_comparison.pdf")
barplot(
  c(auc_combined, auc_individual),
  names.arg = c("Combined", names(auc_individual)),
  las = 2,
  col = c("blue", rep("gray", length(auc_individual))),
  main = "Time-Dependent AUC Comparison",
  ylab = paste("AUC at time =", round(time_point, 2))
)
abline(h = auc_combined, col = "blue", lty = 2)
dev.off()

# --- Save Results ---
write.csv(results_table, "predictive_comparison_results.csv", row.names = FALSE)

# --- Notes for Manuscript ---
cat("
Manuscript Updates:
- Added Methods section describing Cox models, time-dependent AUC, and DeLong’s test.
- Results include Table X with AUC, concordance, and DeLong p-values (see predictive_comparison_results.csv).
- Figure Y shows AUC comparison (see auc_comparison.pdf).
- Combined model AUC:", auc_combined, "vs. individual models (range:", 
min(auc_individual), "-", max(auc_individual), ").
")