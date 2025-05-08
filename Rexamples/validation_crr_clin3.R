#---- bas_cumhazdf_crr_clin3 object created by  call to `crr_fg_train(clin_vars3)`

# --- External Validation on `npxdata_all` ---
library(cmprsk)
library(survival)
library(dplyr)
library(riskRegression)

# --- External Validation on npxdata_all ---

# Extract results from crr_fg_train
fg_fit <- bas_cumhazdf_crr_clin3$fg_fit
bas_cumhazdf_crr <- bas_cumhazdf_crr_clin3$bas_cumhazdf_crr
scaling_params <- bas_cumhazdf_crr_clin3$scaling_params

# Verify bas_cumhazdf_crr
summary(bas_cumhazdf_crr$hazard)  # Expected: ~0.002629076 to 1.442347824
summary(bas_cumhazdf_crr$cif)     # Expected: ~0.002625623 to 0.815776462

# Prepare dataset B (npxdata_all)
X_B <- as.matrix(npxdata_all[, clin_vars3])
X_B_scaled <- scale(X_B, 
                    center = scaling_params$center, 
                    scale = scaling_params$scale)

# Verify data compatibility
if (!all(clin_vars3 %in% colnames(npxdata_all))) {
  stop("clin_vars3 not found in npxdata_all")
}
if (!all(c("event_time", "event_type") %in% colnames(npxdata_all))) {
  stop("event_time or event_type not found in npxdata_all")
}
cat("Event type distribution in npxdata_all:\n")
print(table(npxdata_all$event_type))  # Should include 0, 1, 2

# Compute linear predictor
lin_pred_B <- X_B_scaled %*% fg_fit$coef
lin_pred_B_centered <- lin_pred_B - mean(lin_pred_B)

# Verify linear predictor
summary(lin_pred_B_centered)  # Should be centered, range ~ -3 to 3
if (any(abs(lin_pred_B_centered) > 5)) {
  warning("Large linear predictors detected in lin_pred_B. Check covariate scaling.")
}

# Compute subdistribution survival probabilities and CIF
H_0 <- bas_cumhazdf_crr$hazard
surv_probs <- matrix(NA, nrow = length(H_0), ncol = length(lin_pred_B))
for (i in 1:length(lin_pred_B)) {
  surv_probs[, i] <- exp(-exp(lin_pred_B_centered[i]) * H_0)
}
surv_probs <- as.data.frame(surv_probs)
rownames(surv_probs) <- bas_cumhazdf_crr$time
colnames(surv_probs) <- rownames(X_B)  # Or use IDs from npxdata_all
cif_probs <- 1 - surv_probs

# Verify
summary(surv_probs)  # Range 0 to 1
summary(cif_probs)   # Range 0 to 1

# Calibration: Observed vs. Predicted CIF
cuminc_B <- cuminc(ftime = npxdata_all$event_time, 
                   fstatus = npxdata_all$event_type, cencode = 0)
plot(cuminc_B$"1 1", col = "blue", main = "Observed vs Predicted CIF (Event 1)", 
     xlab = "Time", ylab = "Cumulative Incidence")
mean_cif_pred <- rowMeans(cif_probs)
lines(bas_cumhazdf_crr$time, mean_cif_pred, col = "red", lty = 2)
# Optional: Plot baseline CIF
lines(bas_cumhazdf_crr$time, bas_cumhazdf_crr$cif, col = "green", lty = 3)

# Calibration slope
npxdata_all$lin_pred <- lin_pred_B_centered
cox_calib <- coxph(Surv(event_time, event_type == 1) ~ lin_pred, 
                   data = npxdata_all)
summary(cox_calib)  # Slope near 1 indicates good calibration

# Discrimination (AUC)
auc <- Score(list("Fine-Gray" = cox_calib), 
             formula = Surv(event_time, event_type == 1) ~ 1, 
             data = npxdata_all, times = c(1, 5, 10))
print(auc$AUC)

# Brier score
brier <- Score(list("Fine-Gray" = cox_calib), 
               formula = Surv(event_time, event_type == 1) ~ 1, 
               data = npxdata_all, times = c(1, 5, 10))
print(brier$Brier)