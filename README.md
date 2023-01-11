# Branches

* `2023-01-init`:  5.Cox_tidy.Rmd before dividing into two (independent) parts pertaining to training and validation 


# Notes 

## Categorical covariates
 
 "... if you're interested in interpreting your model or discussing which factors are important after the fact, 
 you're in a weird spot."

## Interactions

https://stackoverflow.com/questions/27580267/how-to-make-all-interactions-before-using-glmnet
https://stats.stackexchange.com/questions/244729/lasso-with-interaction-terms-is-it-okay-if-main-effects-are-shrunk-to-zero


```
library(glmnet)
# Sample data
set.seed(123)
data <- data.frame(matrix(rnorm(9 * 10), ncol = 9))
names(data) <- c(paste0("x", 1:8), "y")
# First step: using .*. for all interactions
f <- as.formula(y ~ .*.)
y <- data$y
# Second step: using model.matrix to take advantage of f
x <- model.matrix(f, data)[, -1]
glmnet(x, y) # x_train, y_train
```

## B-splines

see https://stats.stackexchange.com/questions/373997/spline-regression-with-many-features-in-r

```
library(splines)
f <- as.formula(y ~ .)
xvars <- model.matrix(f, data)[, -1]
colnames(xvars)

# B-spline for x1 
s1 <- bs(data$x1, df = 3)
colnames(s1) <- paste0("s1_", colnames(s1))

# B-spline for x2
s2 <- bs(data$x2, df=4)
colnames(s2) <- paste0("s2_", colnames(s2))

# X matrix for glmnet
xdf <- as.data.frame(cbind(xvars, s1, s2))
xdf$x1 <- NULL  
xdf$x2 <- NULL
x <- as.matrix(xdf)

colnames(x)
```



## The plsmselect package

https://cran.r-project.org/web/packages/plsmselect/vignettes/plsmselect.html

## Links

[Analysis of High-Dimensional Data by Stadler](https://bookdown.org/staedler_n/highdimstats/)
