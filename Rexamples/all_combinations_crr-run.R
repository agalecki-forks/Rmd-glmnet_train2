#-- source("all_combinations_crr-run.R")
library(cmprsk)
library(gtools)
source("../Rfun/all_combinations_crr.R")

###
library(cmprsk)
set.seed(123)

# Create synthetic train_data_all
n <- 156
train_data_all <- data.frame(
  time = rexp(n, rate = 0.1),
  event_type = sample(c(0, 1, 2), n, replace = TRUE, prob = c(94/156, 46/156, 16/156)),
  matrix(rnorm(n * 24), nrow = n, ncol = 24, dimnames = list(NULL, paste0("x", 1:24)))
)

# Introduce signal (mimicking your selected x1, x11, x8, x6)
train_data_all$x1 <- train_data_all$x1 + 0.5 * (train_data_all$event_type == 1)
train_data_all$x11 <- train_data_all$x11 + 0.3 * (train_data_all$event_type == 1)
train_data_all$x8 <- train_data_all$x8 + 0.2 * (train_data_all$event_type == 1)
train_data_all$x6 <- train_data_all$x6 + 0.2 * (train_data_all$event_type == 1)

# Data checks
cat("Event counts:\n")
print(table(train_data_all$event_type))
cat("Number of observations:", nrow(train_data_all), "\n")
cat("Missing values:\n")
print(colSums(is.na(train_data_all)))
cat("Non-positive times:\n")
print(sum(train_data_all$time <= 0))
cor_matrix <- cor(train_data_all[, paste0("x", 1:24)], use = "pairwise.complete.obs")
cat("Max correlation:", max(abs(cor_matrix[upper.tri(cor_matrix)])), "\n")
###
library(gtools)

# Define covariates and xin
covariates <- paste0("x", 2:24)
xin <- c("x1")

# Run all combinations
result <- all_combinations_crr(
  time = "time",
  status = "event_type",
  covariates = covariates,
  xin = xin,
  data = train_data_all,
  criterion = "aic",
  addedx = 3,
  best = 10
)

# Print results
cat("\nTop models summary:\n")
for (i in 1:length(result$top_models)) {
  cat("Model", i, ": Variables =", paste(result$top_models[[i]]$variables, collapse = ", "), 
      ", AIC =", result$top_models[[i]]$criterion, "\n")
}
if (length(result$non_converged) > 0) {
  cat("Non-converged models:", paste(result$non_converged, collapse = "; "), "\n")
}

# Predict CIF for the best model
best_model <- result$top_models[[1]]$model
best_vars <- result$top_models[[1]]$variables
newdata <- data.frame(matrix(0, nrow = 2, ncol = length(best_vars)))
colnames(newdata) <- best_vars
if (length(best_vars) > 0) {
  newdata[2, best_vars[1]] <- 1
}
cov1 <- as.matrix(newdata)
pred <- predict(best_model, cov1 = cov1, times = seq(0, max(train_data_all$time), by = 0.1))
plot(pred, xlab = "Time", ylab = "Cumulative Incidence", main = "CIF for Best Model", col = 1:2)