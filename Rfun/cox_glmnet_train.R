
#' Fit a penalized Cox model using cv.glmnet and extract baseline hazard
#' @param predictors Character vector of predictor names
#' @param event_vars Character vector of event variables (default: c("time", "status"))
#' @param data Data frame containing predictors, event variables, and optional INDEX column
#' @param penalty_factor Numeric vector of penalty factors for predictors
#' @param alpha Elastic net mixing parameter (default: 0.25)
#' @param verbose Logical, print diagnostic outputs (default: TRUE)
#' @param seed Numeric, seed for reproducibility (default: 1234)
#' @return List with cvfit_scaled (cv.glmnet model), bas_cumhazdf (time, hazard),
#'         lin_pred (named numeric vector with INDEX or row names), and scaling_params
cox_glmnet_train <- function(predictors, 
                             event_vars = c("time", "status"), 
                             data = npxdata_all,
                             foldid = NULL, 
                             penalty_factor = NULL, 
                             alpha = 0.25, 
                             verbose = TRUE,
                             seed = 1234) {
  
  # Set seed if provided
  if (!is.null(seed)) set.seed(seed)
  
  # Validate inputs
  if (!all(c(event_vars, predictors) %in% colnames(data))) {
    stop("Some predictors or event variables not found in data: ", 
         paste(setdiff(c(event_vars, predictors), colnames(data)), collapse = ", "))
  }
  if (!all(data[[event_vars[2]]] %in% c(0, 1))) {
    warning("Status should include 0 (censored) and 1 (event)")
  }
  if (is.null(penalty_factor)) {
    penalty_factor <- rep(1, length(predictors))
  } else if (length(penalty_factor) != length(predictors)) {
    stop("penalty_factor length must match number of predictors")
  }
  
  # Define variables
  all_vars <- c(event_vars, predictors)
  if ("INDEX" %in% colnames(data)) {
    all_vars <- c("INDEX", all_vars)
  }
  
  # Subset data, preserving INDEX if present
  data_surv <- data[, all_vars, drop = FALSE]
  
  # Get identifiers (INDEX or row names)
  if ("INDEX" %in% colnames(data_surv)) {
    ids <- data_surv$INDEX
  } else {
    ids <- rownames(data_surv)
    if (is.null(ids) || all(ids == as.character(1:nrow(data_surv)))) {
      warning("No INDEX column or meaningful row names; using row indices as identifiers")
      ids <- as.character(1:nrow(data_surv))
    }
  }
  
  # Standardize covariates
  Xmtx <- as.matrix(data_surv[, predictors, drop = FALSE])
  Xmtx_scaled <- scale(Xmtx)
  if (verbose) {
    cat("Structure of Xmtx_scaled:\n")
    print(str(Xmtx_scaled))
  }
  
  data_surv_scaled <- data.frame(Xmtx_scaled, 
                                 time = data_surv[[event_vars[1]]], 
                                 status = data_surv[[event_vars[2]]])
  if (verbose) {
    cat("Structure of data_surv_scaled:\n")
    print(str(data_surv_scaled))
    cat("Event status distribution:\n")
    print(table(data_surv[[event_vars[2]]]))
  }
  
  # Fit penalized Cox model
  cvfit_scaled <- cv.glmnet(x = Xmtx_scaled, 
                            y = Surv(data_surv_scaled[[event_vars[1]]], 
                                     data_surv_scaled[[event_vars[2]]]),
                            family = "cox", 
                            alpha = alpha, 
                            foldid = foldid,
                            penalty.factor = penalty_factor)
  
  # Compute linear predictors with proper naming
  lin_pred <- predict(cvfit_scaled, newx = Xmtx_scaled, s = "lambda.min", type = "link")
  lin_pred <- as.numeric(lin_pred)  # Convert matrix to numeric vector
  names(lin_pred) <- ids  # Assign INDEX or row names
  if (verbose) {
    cat("Summary of lin_pred:\n")
    print(summary(lin_pred))
    cat("First few names of lin_pred:\n")
    print(head(names(lin_pred)))
    cat("Number of unique lin_pred names:", length(unique(names(lin_pred))), "\n")
  }
  
  # Fit Cox model with offset
  coxph_formula <- as.formula(paste("Surv(", event_vars[1], ",", event_vars[2], ") ~ offset(lin_pred)"))
  coxph_fit1 <- coxph(coxph_formula, data = data.frame(data_surv_scaled, lin_pred = lin_pred))
  
  # Extract cumulative baseline hazard
  bas_cumhazdf <- basehaz(coxph_fit1, centered = FALSE)
  
  # Fix zero hazards and cap outliers
  if (any(bas_cumhazdf$hazard == 0)) {
    if (verbose) message("Zero hazards replaced with minimum positive hazard")
    min_haz <- min(bas_cumhazdf$hazard[bas_cumhazdf$hazard > 0])
    bas_cumhazdf$hazard[bas_cumhazdf$hazard == 0] <- min_haz
  }
  
  # Cap outliers in hazard
  bas_cumhazdf$hazard <- pmin(bas_cumhazdf$hazard, quantile(bas_cumhazdf$hazard, 0.99))
  if (verbose) {
    cat("Summary of bas_cumhazdf$hazard:\n")
    print(summary(bas_cumhazdf$hazard))
  }
  
  # Return results
  return(list(
    cvfit_scaled = cvfit_scaled,
    bas_cumhazdf = bas_cumhazdf,
    lin_pred = lin_pred,
    scaling_params = list(center = attr(Xmtx_scaled, "scaled:center"), 
                          scale = attr(Xmtx_scaled, "scaled:scale"))
  ))
}

# Example usage
#surv_vars <- c("time", "status")
#cvfit_predictors <- c(prot_npx, clin_vars3)  # 21 proteomic + 3 clinical
#table(npxdata_all$status)  # Check distribution of 0, 1
#cox_glmnet_result <- cox_glmnet_train(predictors = cvfit_predictors, 
#                                      penalty_factor = pen, 
#                                      verbose = TRUE)

