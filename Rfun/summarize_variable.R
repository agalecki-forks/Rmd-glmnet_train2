# Function to compute summary statistics for a single numeric variable
summarize_variable <- function(x, var_name) {
  if (!is.numeric(x)) stop("Input 'x' must be numeric")
  if (all(is.na(x))) {
    return(data.frame(
      Variable = var_name,
      n = 0,
      nmiss = length(x),
      mean = NA,
      SD = NA,
      min = NA,
      pct25 = NA,
      pct50 = NA,
      pct75 = NA,
      max = NA,
      stringsAsFactors = FALSE
    ))
  }
  data.frame(
    Variable = var_name,
    n = sum(!is.na(x)),
    nmiss = sum(is.na(x)),
    mean = mean(x, na.rm = TRUE),
    SD = sd(x, na.rm = TRUE),
    min = min(x, na.rm = TRUE),
    pct25 = quantile(x, 0.25, na.rm = TRUE),
    pct50 = quantile(x, 0.50, na.rm = TRUE),
    pct75 = quantile(x, 0.75, na.rm = TRUE),
    max = max(x, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}
