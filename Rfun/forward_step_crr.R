forward_step_crr <- function(time, status, covariates, xin = character(0), data, criterion = "aic", pval_threshold = 0.05) {
  # Inputs:
  #   time: name of the time-to-event variable
  #   status: name of the event status variable (0 = censored, 1 = event of interest, 2 = competing event)
  #   covariates: character vector of candidate covariate names (excluding xin)
  #   xin: character vector of covariates to include by default
  #   data: data frame containing the variables
  #   criterion: "aic" or "pval" for model selection criterion
  #   pval_threshold: threshold for p-value criterion
  
  # Input validation
  if (!time %in% names(data)) stop("Time variable not found in data")
  if (!status %in% names(data)) stop("Status variable not found in data")
  if (!all(covariates %in% names(data))) stop("Some covariates not found in data")
  if (!all(xin %in% names(data))) stop("Some xin covariates not found in data")
  if (anyNA(data[, c(time, status, covariates, xin)])) stop("Missing values detected; please remove or impute")
  if (any(data[[time]] <= 0)) stop("Non-positive time values detected")
  if (!all(data[[status]] %in% c(0, 1, 2))) stop("Status must be 0, 1, or 2")
  if (sum(data[[status]] == 1) < 5) warning("Fewer than 5 events of interest; model may be unstable")
  
  # Initialize model
  current_vars <- xin
  remaining_vars <- setdiff(covariates, xin)
  best_model <- NULL
  best_criterion <- Inf
  
  # Function to compute AIC
  get_aic <- function(model) {
     loglik <- model$loglik  # Use lowercase loglik
    if (is.null(loglik) || is.na(loglik)) {
      warning("Log-likelihood not available; returning Inf")
      return(Inf)
    }
    n_param <- length(model$coef)
    aic <- -2 * loglik + 2 * n_param
    if (is.na(aic) || is.infinite(aic)) {
      warning("Invalid AIC; returning Inf")
      return(Inf)
    }
    return(aic)
  }
  
  cat("Starting forward selection with xin =", xin, "\n")
  cat("Event counts:", table(data[[status]]), "\n")
  
  # Fit initial model with xin (if any)
  if (length(xin) > 0) {
    formula_str <- paste("~", paste(xin, collapse = " + "))
    cat("Initial formula:", formula_str, "\n")
    X <- as.matrix(data[, xin, drop = FALSE])
    best_model <- try(crr(
      ftime = data[[time]],
      fstatus = data[[status]],
      cov1 = X,
      failcode = 1,
      cencode = 0
    ), silent = TRUE)
    if (inherits(best_model, "try-error")) {
      stop("Initial model with xin failed; check predictors for collinearity or variance")
    }
    best_criterion <- get_aic(best_model)
    if (is.infinite(best_criterion)) {
      stop("AIC calculation failed for initial model")
    }
    cat("Initial AIC:", best_criterion, "\n")
  }
  
  while (length(remaining_vars) > 0) {
    cat("Remaining variables:", remaining_vars, "\n")
    best_candidate <- NULL
    best_candidate_criterion <- Inf
    best_candidate_model <- NULL
    
    for (var in remaining_vars) {
      cat("Trying variable:", var, "\n")
      trial_vars <- c(current_vars, var)
      formula_str <- paste("~", paste(trial_vars, collapse = " + "))
      cat("Formula:", formula_str, "\n")
      
      X <- as.matrix(data[, trial_vars, drop = FALSE])
      model <- try(crr(
        ftime = data[[time]],
        fstatus = data[[status]],
        cov1 = X,
        failcode = 1,
        cencode = 0
      ), silent = TRUE)
      
      if (inherits(model, "try-error")) {
        cat("Model failed for", var, "\n")
        next
      }
      
      if (criterion == "aic") {
        crit_value <- get_aic(model)
        if (is.infinite(crit_value)) {
          cat("AIC calculation failed for", var, "\n")
          next
        }
        cat("AIC:", crit_value, "\n")
      } else if (criterion == "pval") {
        model_summary <- summary(model)
        pvals <- model_summary$coef[, "p-value"]
        crit_value <- min(pvals, na.rm = TRUE)
        if (is.na(crit_value) || is.infinite(crit_value)) {
          cat("P-value calculation failed for", var, "\n")
          next
        }
        cat("P-value:", crit_value, "\n")
      } else {
        stop("Invalid criterion. Use 'aic' or 'pval'.")
      }
      
      if (length(crit_value) == 0 || is.na(crit_value)) {
        cat("Invalid criterion value for", var, "\n")
        next
      }
      
      if (crit_value < best_candidate_criterion) {
        best_candidate <- var
        best_candidate_criterion <- crit_value
        best_candidate_model <- model
      }
    }
    
    if (!is.null(best_candidate)) {
      cat("Best candidate:", best_candidate, "with criterion value:", best_candidate_criterion, "\n")
      if (criterion == "aic" && best_candidate_criterion < best_criterion) {
        current_vars <- c(current_vars, best_candidate)
        remaining_vars <- setdiff(remaining_vars, best_candidate)
        best_model <- best_candidate_model
        best_criterion <- best_candidate_criterion
        cat("Added", best_candidate, "to model\n")
      } else if (criterion == "pval" && best_candidate_criterion < pval_threshold) {
        current_vars <- c(current_vars, best_candidate)
        remaining_vars <- setdiff(remaining_vars, best_candidate)
        best_model <- best_candidate_model
        best_criterion <- best_candidate_criterion
        cat("Added", best_candidate, "to model\n")
      } else {
        cat("No improvement, stopping\n")
        break
      }
    } else {
      cat("No valid model found, stopping\n")
      break
    }
  }
  
  return(list(
    selected_model = best_model,
    selected_variables = current_vars,
    criterion_value = best_criterion
  ))
}
