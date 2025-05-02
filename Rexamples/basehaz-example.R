library(survival)
# Example data: time, status, and lp (linear predictor)
data <- data.frame(
  time = c(1, 2, 3, 4, 5),
  status = c(1, 0, 1, 1, 0),
  lp = c(0.1, -0.2, 0.3, 0.0, -0.1)  # Precomputed X*beta
)
# Fit Cox model with offset
fit <- coxph(Surv(time, status) ~ 1 + offset(lp), data = data)


####

# Extract baseline cumulative hazard
base_haz <- basehaz(fit, centered = FALSE)
# Baseline hazard increments
h0 <- diff(c(0, base_haz$hazard))  # Hazard at event times
event_times <- base_haz$time
# Output
print(data.frame(time = event_times, baseline_hazard = h0))