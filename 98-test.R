
rm(list = ls())
f1 <- "./tests/T1.Rdata_cy"
load(file =f1)
ls()
names(T.1)

BICAICglm<- function(fit){
 #-- Based on https://stackoverflow.com/questions/40920051/r-getting-aic-bic-likelihood-from-glmnet
   tmp <- list(fit =fit)
   assign("T.4", tmp, envir =.GlobalEnv)

   #tLL <- fit$null.deviance - deviance(fit)  
   #tLL <- fit$null.deviance - deviance(fit)  
   tLL <- -deviance(fit) # 2*log-likelihood
   ## k <- dim(model.matrix(fit))[2]
   k <- fit$df 
   n <- nobs(fit)
   AICc <- -tLL+2*k+2*k*(k+1)/(n-k-1)
   AIC_ <- -tLL+2*k
   BIC  <-log(n)*k - tLL
   res=cbind(AIC_, BIC, AICc)
   colnames(res)=c("AIC", "BIC", "AICc")
   return(res)
 }
 
mytidy_Surv <- function(cvfit, xnew, ySurv){
  # cv.glmnet fitfamily = cox
   tmp <- list(cvfit= cvfit, xnew=xnew, ySurv= ySurv)
   assign("T.1", tmp, envir = .GlobalEnv)
   mygl <- myglance(cvfit)
   nlmbda <- mygl %>% select(n_lambda) %>% pull()
   idx_min <- mygl[, "index_min"] %>% pull()
   idx_1se <- mygl[, "index_1se"] %>% pull()
   cvtd <- mytidy(cvfit) 
   # print(cvtd)
   mincv <- rep("-", nlmbda)
   mincv[idx_1se:idx_min] <- "+"
   mincv[idx_min] <- "min>"
   mincv[idx_1se] <- "<1se"
   tmp <- list(mygl= mygl, cvtd = cvtd, idx = c(idx_1se, idx_min), mincv=mincv)
   assign("T.2", tmp, envir = .GlobalEnv)

   # print(length(mincv))
   cvtd$mincv <- mincv
 
   gfit <- cvfit$glmnet.fit
 
   td_fit <- mytidy(gfit) %>% select(-c(step, lambda)) # with nested beta
   pred     <- predict(gfit, newx = xnew)
   tmp <- list(pred = pred, gfit = gfit)
   assign("T.3", tmp, envir = .GlobalEnv)

   C_index  <-  Cindex(pred, ySurv)
 
   info     <- BICAICglm(gfit)
   info_tbl <- as_tibble(info)
   td_fit$Cindex <- C_index
   bind_cols(cvtd, td_fit, info_tbl)
}

 if (!require('pacman')) install.packages('pacman', repos = "http://cran.us.r-project.org")
 library(pacman)
 
 pacman::p_load(
   readxl,       # load excel file
   glmnet,       # Lasso and Elastic-Net Regularized Generalized Linear Models
   rmdformats,   # rmd formats
   rmarkdown,    # rmarkdown
   here,         # File locator
   skimr,        # get overview of data
   tidyverse,    # data management + ggplot2 graphics 
   tidymodels,   #for modeling and machine learning using tidyverse principles
   janitor,      # adding totals and percents to tables
   flextable,    # format the output table
   gtsummary,    # calculate summary statistics and format the results
   sjPlot,       # correlation matrix
   purrr,        # enhances R functional programming (FP) toolki 
   tidyr,        #Tools to help to create tidy data
   ggplot2,      #Plot the results
   glmnetUtils,  #Glmnet models for multiple alpha
   coefplot,     # Plotting Model Coefficients
   survival,     #survival model 
   survivalROC,  #survivalROC
   writexl       #write excel
   )

library(utilsag)
# session Info

sessionInfo()

# unpack T.1 

names(T.1)
cvfit1 <- T.1[["cvfit"]]
x  <- T.1[["xnew"]]
ySurv <- T.1[["ySurv"]]
td_cvfit1 <- mytidy_Surv(cvfit1, x, ySurv)
td_cvfit1
ls()


