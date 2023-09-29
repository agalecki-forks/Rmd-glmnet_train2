# source("18.Cox-multiple-alpha-run.R")
# Clears Global environment
 rm(list=ls())
 T.1 <- T.2 <- T.3 <- NULL  # Use 99_test.R to examine these objects
 
 # Auxiliary functions
  
 names3_aux <- function(basenm, postfix){ 
    nmRmd <- paste0(basenm, ".Rmd")
    nmR <- paste0("./purl/", basenm, ".Rprog")
    nm_out <- paste0( "./out/18out/html/", basenm, postfix)
    c(nmRmd = nmRmd, nmR = nmR, nm_out = nm_out)  
 }
 
 BICAICglm<- function(fit){
 #-- Based on https://stackoverflow.com/questions/40920051/r-getting-aic-bic-likelihood-from-glmnet
   tmp <- list(fit =fit)
   # assign("T.4", tmp, envir =.GlobalEnv)

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
  # cv.glmnet fit, family = cox
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
 
   fit <- cvfit$glmnet.fit
 
   td_fit <- mytidy(fit) %>% select(-c(step, lambda)) # with nested beta
   pred     <- predict(fit, newx = xnew)
   tmp <- list(pred = pred, fit = fit)
   assign("T.3", tmp, envir = .GlobalEnv)

   ##-C_index  <-  Cindex(pred, ySurv) # corrected
   C_index <- apply(pred, 2, Cindex, y=ySurv)
   info     <- BICAICglm(fit)
   info_tbl <- as_tibble(info)
   td_fit$Cindex <- C_index
   bind_cols(cvtd, td_fit, info_tbl)
}
# mytidy_Surv(cvfit0, x, ySurv)


 
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
   tidyr,        # Tools to help to create tidy data
   ggplot2,      # Plot the results
   glmnetUtils,  # Glmnet models for multiple alpha
   coefplot,     # Plotting Model Coefficients
   survival,     # survival model 
   survivalROC,  # survivalROC
   writexl       # write excel
   )

# session Info
session_Info <- sessionInfo()
alphax_vals <- c(0, 0.25, 0.5, 0.75, 1)


for (alphai in 1:length(alphax_vals)){

 alphax <- alphax_vals[alphai]

Rdata_out <- paste0("./out/18out/a_", alphax, "/18-cox-multiple-alpha.Rdata")

save(session_Info, file = Rdata_out)

# Set value of alphax !!!!


# process -init.Rmd
bnm <- "18.Cox-multiple-alpha"   # Basename of Rmd file (do not change it)
bnm_init <- paste0(bnm, "-init.Rmd")
out_init <- paste0("./out/18out/a_", alphax,"/", bnm, "-init")
print(out_init)
parms <- list(alphax = alphax)

rmarkdown::render(bnm_init, "all", output_file = out_init, params = parms)
# saved objects for later use: cvfit1, td_cvfit1 
# saved tibbles for later use: td_cv1, beta1

# --- Loop over columns of x matrix

istart  <- length(clin_vars) + 1
iend    <- length(clin_vars) + length(prot_npx)
# iend    <- length(clin_vars) +2  # !!! atg for testing 
loopii <- istart:iend
## loopii <- c(9, 23) # !! atg for testing
len <- length(loopii)

all_cvfit  <- vector("list", length    = len )
all_betas  <- vector("list", length    = len )      # multiple rows
all_cvtidy    <- vector("list", length = len )   
all_cvtd      <- vector("list", length = len )   


for (ix in 1:len){
 ii <- loopii[ix]
 ## ac <- gsub("[.]", "_", paste0("-", a)) # . -> _
 message ("--- alphax=", alphax, ", ii=", ii,": Rmd for covariate= ", xvars[ii], " processed")
 paramsi <- list(alphax = alphax, ei = ii)   # index of covariate excluded
 nmsj <- names3_aux(bnm, ii)
 knitr::purl(nmsj["nmRmd"], output = nmsj["nmR"])
 out_nm <- paste0("./out/18out/a_", alphax,"/html/", bnm, ii)
 rmarkdown::render(nmsj["nmRmd"], "all", output_file = out_nm , params = paramsi)
 # print(cvfiti)
 # save objects : cvfit, td_cvfit 
 # save tibbles: td_cv, betai

 #-  cvfit, td_cvfit, betai
 all_cvfit[[ix]]  <- cvfit
 all_betas[[ix]]  <- betai   
 all_cvtidy[[ix]]  <- td_cvfit
 all_cvtd[[ix]]   <- td_cv

}

x_nms0 <- paste("x", loopii,"_", colnames(x)[loopii], sep ="")
x_nms  <- c("ALLx", x_nms0) 
res_cvfit <- append(all_cvfit,  list(cvfit1),    after = 0)
res_betas <- append(all_betas,  list(beta1),     after = 0)
res_cvtidy <-append(all_cvtidy, list(td_cvfit1), after = 0)
res_cvtd   <-append(all_cvtd,   list(td_cv1),    after = 0) # Not saved

names(res_cvfit) <- x_nms
names(res_betas) <- x_nms
names(res_cvtidy) <- x_nms
names(res_cvtd) <- x_nms
save(alphax, res_cvfit, res_betas, res_cvtidy, file = Rdata_out)
# load(Rdata_out,verbose = TRUE)

#-- Create xlsx

xlsx_path <- paste0("./out/18out/a_", alphax, "/res_betas.xlsx")
write_xlsx(res_betas, xlsx_path)
xlsx_path <- paste0("./out/18out/a_", alphax, "/res_cvtidy.xlsx")
write_xlsx(res_cvtidy, xlsx_path)
} # for alphai 
