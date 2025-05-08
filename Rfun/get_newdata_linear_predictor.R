get_newdata_linear_predictor <- function(cvfit, newdata) {
  if (!inherits(cvfit, "cv.glmnet")) {
    stop("cvfit must be a cv.glmnet object")
  }
  
  # Check if newdata has the correct number of columns
  expected_cols <- nrow(cvfit$glmnet.fit$beta)
  if (ncol(newdata) != expected_cols) {
    stop("newdata must have ", expected_cols, " columns, but has ", ncol(newdata))
  }
  
  # Ensure newdata is a matrix
  newdata_matrix <- as.matrix(newdata)
  if (any(is.na(newdata_matrix))) {
    stop("newdata contains missing values")
  }
  
  # Predict linear predictor using lambda.min
  lp <- predict(cvfit, newx = newdata_matrix, s = "lambda.min", type = "link")
  
  return(as.vector(lp))
}