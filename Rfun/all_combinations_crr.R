
all_combinations_crr <- function(time, status, covariates, xin = character(0), data, criterion = "aic", pval_threshold = 0.05, addedx = 10, best = 10) {
  # Inputs:
  #   time: name of the time-to-event variable
  #   status: name of the event status variable (0 = censored, 1 = event of interest, 2 = competing event)
  #   covariates: character vector of candidate covariate names (excluding xin)
  #   xin: character vector of covariates to include by default
  #   data: data frame containing the variables
  #   criterion: "aic" or "pval" for model selection criterion
  #   pval_threshold: threshold for p-value criterion (used if criterion = "pval")
  #   addedx: maximum number of covariates to add from covariates (default 10)
  #   best: number of top models to return (default 10)
  
  # Input validation
  if (!time %in% names(data)) stop("Time variable not found in data")
  if (!status %in% names(data)) stop("Status variable not found in data")
  if (!all(covariates %in% names(data))) stop("Some covariates not found in data")
  if (!all(xin %in% names(data))) stop("Some xin covariates not found in data")
  if (anyNA(data[, c(time, status, covariates, xin)])) stop("Missing values detected; please remove or impute")
  if (any(data[[time]] <= 0)) stop("Non-positive time values detected")
  if (!all(data[[status]] %in% c(0, 1, 2))) stop("Status must be 0, 1, or 2")
  if (sum(data[[status]] == 1) < 5) warning("Fewer than 5 events of interest; model may be unstable")
  if (addedx < 1 || addedx > length(covariates)) stop("addedx must be between 1 and length(covariates)")
  if (best < 1) stop("best must be at least 1")
  
  # Initialize
  model_results <- list()
  non_converged <- character()
  
  # Function to compute AIC
  get_aic <- function(model) {
    loglik <- model$loglik
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
  
  cat("Fitting models with xin =", xin, "and up to", addedx, "added covariates\n")
  cat("Event counts:", table(data[[status]]), "\n")
  
  # Fit initial model with xin (if any)
  if (length(xin) > 0) {
    formula_str <- paste("~", paste(xin, collapse = " + "))
    cat("Initial formula:", formula_str, "\n")
    X <- as.matrix(data[, xin, drop = FALSE])
    model <- try(crr(
      ftime = data[[time]],
      fstatus = data[[status]],
      cov1 = X,
      failcode = 1,
      cencode = 0
    ), silent = TRUE)
    
    if (inherits(model, "try-error") || !model$converged) {
      stop("Initial model with xin failed or did not converge; check predictors for collinearity or variance")
    }
    
    crit_value <- if (criterion == "aic") {
      get_aic(model)
    } else {
      pvals <- summary(model)$coef[, "p-value"]
      min(pvals, na.rm = TRUE)
    }
    
    if (is.infinite(crit_value) || is.na(crit_value)) {
      stop("Criterion calculation failed for initial model")
    }
    
    model_results[[1]] <- list(
      variables = xin,
      model = model,
      criterion = crit_value
    )
    cat("Initial", criterion, ":", crit_value, "\n")
  }
  
  # Generate all possible combinations of 1 to addedx covariates
  for (k in 1:min(addedx, length(covariates))) {
    cat("Trying combinations of", k, "additional covariates\n")
    combos <- combinations(n = length(covariates), r = k, v = covariates)
    
    for (i in 1:nrow(combos)) {
      trial_vars <- c(xin, combos[i, ])
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
      
      if (inherits(model, "try-error") || !model$converged) {
        non_converged <- c(non_converged, paste(combos[i, ], collapse = ", "))
        cat("Model failed or did not converge for", paste(combos[i, ], collapse = ", "), "\n")
        next
      }
      
      if (criterion == "aic") {
        crit_value <- get_aic(model)
        if (is.infinite(crit_value)) {
          cat("AIC calculation failed for", paste(combos[i, ], collapse = ", "), "\n")
          next
        }
        cat("AIC:", crit_value, "\n")
      } else if (criterion == "pval") {
        model_summary <- summary(model)
        pvals <- model_summary$coef[, "p-value"]
        crit_value <- min(pvals, na.rm = TRUE)
        if (is.na(crit_value) || is.infinite(crit_value)) {
          cat("P-value calculation failed for", paste(combos[i, ], collapse = ", "), "\n")
          next
        }
        cat("P-value:", crit_value, "\n")
      } else {
        stop("Invalid criterion. Use 'aic' or 'pval'.")
      }
      
      if (length(crit_value) == 0 || is.na(crit_value)) {
        cat("Invalid criterion value for", paste(combos[i, ], collapse = ", "), "\n")
        next
      }
      
      model_results[[length(model_results) + 1]] <- list(
        variables = trial_vars,
        model = model,
        criterion = crit_value
      )
    }
  }
  
  # Sort and select top 'best' models
  if (length(model_results) == 0) {
    stop("No valid models found")
  }
  
  crit_values <- sapply(model_results, function(x) x$criterion)
  top_indices <- order(crit_values)[1:min(best, length(model_results))]
  top_models <- model_results[top_indices]
  
  # Summarize results
  cat("\nTop", length(top_models), "models (sorted by", criterion, "):\n")
  for (i in 1:length(top_models)) {
    cat("Model", i, ": Variables =", paste(top_models[[i]]$variables, collapse = ", "), 
        ",", criterion, "=", top_models[[i]]$criterion, "\n")
  }
  
  if (length(non_converged) > 0) {
    cat("\nNon-converged models:", length(non_converged), "\n")
    cat(paste(non_converged, collapse = "; "), "\n")
  }
  
  return(list(
    top_models = top_models,
    non_converged = non_converged
  ))
}
