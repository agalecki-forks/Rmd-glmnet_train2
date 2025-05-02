#-- source("forward_step_crr-run.R")
library(cmprsk)
source("../Rfun/forward_step_crr.R")

#####
library(timereg)
set.seed(123)

# Create synthetic train_data_all
n <- 156
train_data_all <- data.frame(
  time = rexp(n, rate = 0.1),  # Exponential times
  event_type = sample(c(0, 1, 2), n, replace = TRUE, prob = c(91/156, 51/156, 14/156)),
  matrix(rnorm(n * 24), nrow = n, ncol = 24, dimnames = list(NULL, paste0("x", 1:24)))
)

# Adjust predictors to introduce some signal
train_data_all$x1 <- train_data_all$x1 + 0.5 * (train_data_all$event_type == 1)  # Stronger effect for x1
train_data_all$x2 <- train_data_all$x2 + 0.2 * (train_data_all$event_type == 1)  # Weaker effect for x2

# Verify event counts
cat("Event counts:\n")
print(table(train_data_all$event_type))
cat("Number of observations:", nrow(train_data_all), "\n")

# Define covariates and xin
covariates <- paste0("x", 1:24)
xin <- c("x1")  # Include x1 by default

# Data checks
cat("Missing values:\n")
print(colSums(is.na(train_data_all)))
cat("Non-positive times:\n")
print(sum(train_data_all$time <= 0))
cor_matrix <- cor(train_data_all[, paste0("x", 1:24)], use = "pairwise.complete.obs")
cat("Max correlation:", max(abs(cor_matrix[upper.tri(cor_matrix)])), "\n")

# Run forward stepwise selection
result <- forward_step_crr(
  time = "time",
  status = "event_type",
  covariates = covariates,
  xin = xin,
  data = train_data_all,
  criterion = "aic"
)

# Print results
cat("Selected variables:", result$selected_variables, "\n")
cat("AIC:", result$criterion_value, "\n")
summary(result$selected_model)

# Predict CIF
newdata <- data.frame(matrix(0, nrow = 2, ncol = length(result$selected_variables)))
colnames(newdata) <- result$selected_variables
if (length(result$selected_variables) > 0) {
  newdata[2, result$selected_variables[1]] <- 1
}
cov1 <- as.matrix(newdata)
pred <- predict(result$selected_model, cov1 = cov1, times = seq(0, max(train_data_all$time), by = 0.1))

plot(pred, xlab = "Time", ylab = "Cumulative Incidence", main = "CIF for Event 1", col = 1:2)

















































