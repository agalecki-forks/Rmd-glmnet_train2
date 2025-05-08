#' Fit a Fine-Gray competing risks model and extract baseline hazard and CIF
#' @param predictors Character vector of predictor names
#' @param event_vars Character vector of event variables (default: c("event_time", "event_type"))
#' @param data Data frame containing predictors and event variables
#' @param verbose Logical, print diagnostic outputs (default: TRUE)
#' @return List with fg_fit (crr model), bas_cumhazdf_crr (time, hazard, cif), 
#'         lin_pred_A (linear predictors for training data), and scaling_params
crr_fg_train <- function(predictors, 
                         event_vars = c("event_time", "event_type"), 
                         data = npxdata_all, 
                         verbose = TRUE) {
  # Validate inputs
  if (!all(c(event_vars, predictors) %in% colnames(data))) {
    stop("Some predictors or event variables not found in data")
  }
  if (!all(c(0, 1, 2) %in% data[[event_vars[2]]])) {
    warning("event_type should include 0 (censored), 1 (event of interest), 2 (competing)")
  }
 
all_vars <- c(event_vars, predictors)

# Subset data
data_surv <- data %>% select(all_of(all_vars))

# Standardize covariates (to prevent large lin_pred)
X_A <- as.matrix(data_surv[, predictors])
X_A_scaled <- scale(X_A)
cat("str X_A_scaled")
print(str(X_A_scaled))

data_surv_scaled <- data.frame(X_A_scaled, event_time = data_surv$event_time, 
                               event_type = data_surv$event_type)
if (verbose){
 cat("Structure of data_surv_scaled")
 print(str(data_surv_scaled))
 cat("Event type distribution:\n")
 print(table(data_surv$event_type))
} 

# Fit Fine-Gray model
fg_fit <- crr(ftime = data_surv_scaled$event_time, 
              fstatus = data_surv_scaled$event_type, 
              cov1 = as.matrix(data_surv_scaled[, predictors]), 
              failcode = 1, cencode = 0)

# Extract cumulative baseline subdistribution hazard
zero_cov <- matrix(0, nrow = nrow(data_surv_scaled), ncol = length(predictors))
colnames(zero_cov) <- predictors

times <- sort(unique(data_surv_scaled$event_time))
pred_baseline <- predict(fg_fit, cov1 = zero_cov, time = times)
if (verbose){
 cat("Structure of pred_baseline \n")
 print(str(pred_baseline))
}

# Create bas_cumhazdf_crr (generalized for matrix or list output)
if (is.matrix(pred_baseline)) {
  if (ncol(pred_baseline) < 2) stop("pred_baseline has insufficient columns")
  
  # Matrix output: first column is times, others are CIF values
  cif_baseline <- rowMeans(pred_baseline[, 2:ncol(pred_baseline)])
 
  bas_cumhazdf <- data.frame(
    time = pred_baseline[, 1],
    hazard = -log(1 - cif_baseline),
    cif  = cif_baseline   # baseline CIF values.
     )
} else {
  # List output: use $time and $cif
  cif_baseline <- rowMeans(pred_baseline$cif)

  bas_cumhazdf <- data.frame(
    time = pred_baseline$time,
    hazard = -log(1 - cif_baseline),
    cif = cif_baseline
  )
}
 if (verbose) {
    cat("Summary of bas_cumhazdf$hazard:\n")
    print(summary(bas_cumhazdf$hazard))
    cat("Summary of bas_cumhazdf$cif:\n")
    print(summary(bas_cumhazdf$cif))
  }
 
  
# Fix zeros and outliers in hazard
  min_haz <- min(bas_cumhazdf$hazard[bas_cumhazdf$hazard > 0])
  if (any(bas_cumhazdf$hazard == 0)) {
    if (verbose) message("Zero hazards replaced with minimum positive hazard")
    bas_cumhazdf$hazard[bas_cumhazdf$hazard == 0] <- min_haz
  }
  bas_cumhazdf$hazard <- pmin(bas_cumhazdf$hazard, quantile(bas_cumhazdf$hazard, 0.99))

# Extract linear predictor for dataset A
  lin_pred_A <- as.matrix(data_surv_scaled[, predictors]) %*% fg_fit$coef
  if (verbose) {
    cat("Summary of lin_pred_A:\n")
    print(summary(lin_pred_A))
  }

 return(list(
  fg_fit = fg_fit,
  bas_cumhazdf_crr = bas_cumhazdf,
  lin_pred = lin_pred_A,
  scaling_params = list(center = attr(X_A_scaled, "scaled:center"), 
                      scale  = attr(X_A_scaled, "scaled:scale"))
 ))  
} # End of function
