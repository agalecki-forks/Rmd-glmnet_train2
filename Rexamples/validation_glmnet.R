# cox_glmnet_result object created by:
# cox_glmnet_result <- cox_glmnet_train(predictors = cvfit_predictors, 
#                                      penalty_factor = pen, 
#                                      verbose =TRUE)

# --- External Validation on `npxdata_all` ---

# Extract results
cvfit_scaled <- cox_glmnet_result$cvfit_scaled
bas_cumhazdf <- cox_glmnet_result$bas_cumhazdf
lin_pred_A <- cox_glmnet_result$lin_pred
scaling_params <- cox_glmnet_result$scaling_params

# Prepare dataset B
X_B <- as.matrix(npxdata_all[, cvfit_predictors])
X_B_scaled <- scale(X_B, 
                    center = scaling_params$center, 
                    scale = scaling_params$scale)

# Compute linear predictor
lin_pred_B <- predict(cvfit_scaled, newx = X_B_scaled, s = "lambda.min", type = "link")
lin_pred_B_centered <- lin_pred_B - mean(lin_pred_B)

# Verify
summary(lin_pred_B_centered)  # Should be centered, range ~ -3 to 3
if (any(abs(lin_pred_B_centered) > 5)) {
  warning("Large linear predictors detected in lin_pred_B. Check covariate scaling.")
}

# Compute survival probabilities
H_0 <- bas_cumhazdf$hazard
surv_probs <- matrix(NA, nrow = length(H_0), ncol = length(lin_pred_B))
for (i in 1:length(lin_pred_B)) {
  surv_probs[, i] <- exp(-exp(lin_pred_B_centered[i]) * H_0)
}
surv_probs <- as.data.frame(surv_probs)
rownames(surv_probs) <- bas_cumhazdf$time
colnames(surv_probs) <- rownames(X_B)  # Or use IDs from npxdata_all

# Verify
summary(surv_probs)  # Range 0 to 1

# Calibration
km_B <- survfit(Surv(time, status) ~ 1, data = npxdata_all)
plot(km_B, col = "blue", main = "Observed vs Predicted Survival", 
     xlab = "Time", ylab = "Survival Probability")
mean_surv_pred <- rowMeans(surv_probs)
lines(bas_cumhazdf$time, mean_surv_pred, col = "red", lty = 2)

npxdata_all$lin_pred <- lin_pred_B_centered
cox_calib <- coxph(Surv(time, status) ~ lin_pred, data = npxdata_all)
summary(cox_calib)  # Slope near 1

# Discrimination
auc <- Score(list("Cox" = cox_calib), 
             formula = Surv(time, status) ~ 1, 
             data = npxdata_all, times = c(1, 5, 10))
print(auc$AUC)

# Brier score
brier <- Score(list("Cox" = cox_calib), 
               formula = Surv(time, status) ~ 1, 
               data = npxdata_all, times = c(1, 5, 10))
print(brier$Brier)