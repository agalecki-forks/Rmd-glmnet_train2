# source("15.Cox-multiple-run.R")
# Clears Global environment
 rm(list=ls())
 Rdata_out <- "15.test.Rdata"

 
 # Auxiliary functions
 
 
 names3_aux <- function(basenm, postfix){ 
    nmRmd <- paste0(basenm, ".Rmd")
    nmR <- paste0("./purl/", basenm, ".Rprog")
    nm_out <- paste0( "./out/15out/html/", basenm, postfix)
    c(nmRmd = nmRmd, nmR = nmR, nm_out = nm_out)  
 }
 
 BICAICglm<- function(fit){
 #-- Based on https://stackoverflow.com/questions/40920051/r-getting-aic-bic-likelihood-from-glmnet
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
   mygl <- myglance(cvfit)
   nlmbda <- mygl %>% select(n_lambda) %>% pull()
   idx_min <- mygl[, "index_min"] %>% pull()
   idx_1se <- mygl[, "index_1se"] %>% pull()
   
   cvtd <- mytidy(cvfit) 
   print(cvtd)
   mincv <- rep(".", nlmbda)
   mincv[idx_min] <- "min"
   mincv[idx_1se] <- "1se"
   print(length(mincv))
   cvtd$mincv <- mincv
   
   fit <- cvfit$glmnet.fit
   td_fit <- mytidy(fit) %>% select(-c(step, lambda)) # with nested beta
   pred     <- predict(fit, newx = xnew)
   C_index  <-  Cindex(pred, ySurv)
   info     <- BICAICglm(fit)
   info_tbl <- as_tibble(info)
   td_fit$Cindex <- C_index
   bind_cols(cvtd, td_fit, info_tbl)
}
# mytidy_Surv(fit0, x, ySurv)


 
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

# session Info
session_Info <- sessionInfo()



save(session_Info, file = Rdata_out)


bnm <- "15.Cox-multiple"   # Basename of Rmd file (do not change it)
bnm_init <- paste0(bnm, "-init.Rmd")
rmarkdown::render(bnm_init, "all")

# --- Loop over columns of x matrix

istart  <- length(clin_vars)+ 1
iend    <- length(clin_vars) + length(prot_npx)
iend    <- length(clin_vars) +2  # !!! atg for testing 
len <- iend - istart + 1

all_cvfit  <- vector("list", length =len )
all_betas  <- vector("list", length =len )   # multiple rows
all_res    <- vector("list", length =len )   # one row




for (ii in istart:iend){
 ## ac <- gsub("[.]", "_", paste0("-", a)) # . -> _
 message ("--- Rmd for covariate= ", xvars[ii], " processed")
 paramsi <- list(ei = ii)   # index of covarite excluded
 nmsj <- names3_aux(bnm, ii)
 knitr::purl(nmsj["nmRmd"], output = nmsj["nmR"])
 rmarkdown::render(nmsj["nmRmd"], "all", output_file = nmsj["nm_out"], params = paramsi)
 # print(cvfiti)
 all_cvfit[[ii]]  <- cvfiti
 all_betas[[ii]]  <- betai    
 all_res[[ii]]    <- resi
}

x_nms <- paste("x_", 1:len)
names(all_cvfit) <- x_nms
save(session_Info, all_cvfit, all_betas, all_res, file = Rdata_out)
