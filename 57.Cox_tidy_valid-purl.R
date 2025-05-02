## ----setup, include=FALSE-------------------------------------------
knitr::opts_chunk$set(echo = TRUE, comment="#>")


## ----data, include = FALSE------------------------------------------
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
  tidymodels,   #for modeling and machine learning using tidyverse principles
  survivalROC   #survivalROC
  )


## ----utilsag-version-2test, include = FALSE-------------------------
uag_ver0 <- "0.2.1" # utilsag version tested to execute this script


## ----install-tested-utilsag-child, child = "_install-tested-utilsag.Rmd", include = FALSE----

## ----utilsag-info, include = FALSE----------------------------------
install.packages("penAFT")
.libPaths()
uag_path <- system.file(package = "utilsag")
uag_ver  <- if (uag_path != ""){
     as.character(packageVersion("utilsag"))} else ""


## ----install-tested-utilsag, include=FALSE--------------------------
uag_ref0 <- paste0("version-", uag_ver0)

if("utilsag" %in% (.packages())) detach("package:utilsag", unload=TRUE)

if (uag_ver != uag_ver0){
   devtools::install_github("agalecki/utilsag", ref = uag_ref0)
}



## ----utilsag-lodaed-------------------------------------------------
library(utilsag)
mod_lbl <-"M3"
alphax <- 0.25
survSplit_cut <- 10



## ----prepare-data, echo=TRUE, message=FALSE, warning=FALSE----------
#read the data that is stored under the `validation_input` folder

# Change the filename below, if needed
dt <- read_excel("./validation_input/data_validation.xlsx")

dt <- dt  %>% rename(time = FU_TIME, status = CASE_CONTROL) %>% 
  mutate(log10_DU_ACR = log10(DU_ACR))  %>% filter(time>0) %>% filter(BL_eGFR >= 45)

dim(dt) # Number of rows and columns in the input data


## ----Variables-used-Cox---------------------------------------------
prot_npx <- c("KIM1.npx","SYND1.npx","IL.1RT1.npx",   "WFDC2.npx", "CD27.npx",
              "TNFRSF10A.npx","LAYN.npx","PVRL4.npx", "EDA2R.npx","TNFRSF4.npx",
              "GFR_alpha_1_npx","TNF.R1.npx","PI3_npx", "EFNA4.npx","TNF.R2.npx" ,
              "DLL1.npx", "TNFRSF6B.npx", "CD160.npx", "EPHA2.npx","RELT.npx",
              "LTBR.npx") 
surv_vars <-  c("time", "status")
clin_vars <- c("BL_eGFR","B_HBA1C_PRC","log10_DU_ACR","SEX","AGE_TL")
clin_vars3 <-  c("BL_eGFR","B_HBA1C_PRC","log10_DU_ACR")
xvars <- c(clin_vars3, prot_npx)


## ----select-mvars, include =FALSE-----------------------------------

if (mod_lbl == "M0") mvars <- c(surv_vars, clin_vars) 
if (mod_lbl == "M1") mvars <- c(surv_vars, prot_npx, clin_vars) 
if (mod_lbl == "M2") mvars <- c(surv_vars, prot_npx, clin_vars) 
if (mod_lbl == "M3") mvars <- c(surv_vars, prot_npx, clin_vars3) 



## ----vars vector----------------------------------------------------
mvars


## ----data-prep1-----------------------------------------------------
data <- dt  %>% select(all_of(mvars)) 
dim(data)         


## ----data-prep2-----------------------------------------------------
# drop any records with NA
test_data <- data %>% drop_na()    

# Total sample size
nrow(test_data)


## ----load-train-objects---------------------------------------------
fpath1 <- paste0("./validation_input/17-cox-0.25", mod_lbl, ".Rdata")
load (file = fpath1, verbose = TRUE)
#- fpath2 <- paste0("./save/5.Cox_tidy2_", mod_lbl, ".Rdata")
# save(mod_selected, cva_pfit, mod_all, sel3long, file = fpath2)



## ----hyperparametrs-selected----------------------------------------
alphax


## ----survSplit------------------------------------------------------
message("survSplit_cut =", survSplit_cut)

temp <- survSplit(Surv(time,  status) ~ ., data = test_data, cut = survSplit_cut,
                  episode = "epsd")
dim(test_data)
test_data15 <- subset(temp, epsd == 1)  # only the first ?? years
dim(test_data15)
test_data_saved <- test_data
test_data    <-  test_data15


## ----x-surv-test----------------------------------------------------
# - dtx_test <- subset(test_data, select=-c(time,status))
dtx_test <- test_data %>% select(all_of(xvars)) 
x_test <- model.matrix(~0 +., data=dtx_test)
dim(x_test)
colnames(x_test)
# Create Surv object
y_test <- data.matrix(test_data[,c("time", "status")])
ySurv_test <- survival::Surv(y_test[,"time"], y_test[, "status"])


## ----auxx-----------------------------------------------------------
cens10 <- ifelse(y_test[,"time"]  < 10 & y_test[, "status"]== 0, 1, 0)

table(cens10)
ttx <- cbind(y_test[,"time"], y_test[, "status"], cens10)
colnames(ttx) <- c("time", 'status', "cens10")
ttx[1:10,]
Cp10_test <- sum(cens10)/length(cens10) 
nnc <- sum(cens10 ==0) # number of non-censored obs
nev <- sum(cens10 ==0 & y_test[, "status"]==1)
prev10_test <- nev/nnc  #prevalence
Cp_test <- Cp10_test
prev_test <- prev10_test

# Cleanup
Cp10_test <- prev10_test <- NULL


## ----cva-glmnet-fit-------------------------------------------------
names(res_cvfit)
cvfit <- res_cvfit$ALLx

pfit_aopt <- cvfit$glmnet
class(pfit_aopt)
names(pfit_aopt)
pfit_aopt$lambda

lmbda <- pfit_aopt$lambda
sel <- c(1, 15, 18, 24, 31) # 30 -> 31
nms <- paste0("step", sel)
lmbda_sel <-lmbda[sel]
names(lmbda_sel) <- nms
lmbda_sel
lmbda_opt <- lmbda_sel[1]
lmbda_opt


## ----calc-pred2-----------------------------------------------------
predM_lpmtx <- predict(pfit_aopt, newx = x_test, type = "link") # Matrix


## ----C-index-Mod----------------------------------------------------
Cindex_lmbda <- apply(predM_lpmtx, 2, Cindex, y = ySurv_test) # Multiple lambda
length(Cindex_lmbda)

Cindex_lmbda[sel] 

#Cindex(predM_lp, ySurv_test)   # For optimal lambda


## ----Cindex1--------------------------------------------------------
cindex_sel <- sapply(sel, FUN = function(sx){
   pred1 <- predM_lpmtx[, sx]
   pred1m <- -pred1
   obj   <- concordance(ySurv_test ~ pred1m)
   cidx  <- obj$concordance
   stder <- sqrt(obj$var)
   res <- tibble(sel = sx, cindex= cidx, stderr = stder)
   return(res)
})

cidx <- t(cindex_sel)
rownames(cidx) <- names(Cindex_lmbda[sel]) 
cidx


## ----res-list-------------------------------------------------------
lmbda_sel
len <- length(lmbda_sel)
survROC_list <- vector(mode ="list", length = len)
names(survROC_list) <- names(lmbda_sel)
str(survROC_list)


## ----score-details--------------------------------------------------
tmp1 <- apply(x_test, 2, mean)
tmp2 <- tmp1
tmp2["BL_eGFR"] <- tmp1["BL_eGFR"]+10
x_toy <- rbind(tmp1,tmp2)
predM_lp <- as.vector(predict(pfit_aopt, newx = x_toy, type = "link", s = lmbda_opt))
res <-cbind(x_toy[, "BL_eGFR"], predM_lp) 
colnames(res) <- c("BL_eGFR", "lp_score")
res


## ----step1-roc, echo = FALSE----------------------------------------
step_sel <- "step1"


## ----surv-roc-------------------------------------------------------
dataM_test <- test_data
dataM_test$predM_lp <- NULL


lmbda_opt <-  lmbda_sel[names(lmbda_sel) == step_sel]
predM_lp <- as.vector(predict(pfit_aopt, newx = x_test, type = "link", s = lmbda_opt))
# predM_lp <- - predM_lp # !!!
predM_lp[1:10]
summary(predM_lp)

                
## Augment `dataM_test` with linear predictor
 dataM_test$predM_lp <- predM_lp

# Evaluate every 2.5 years
 
tx <- 2.5* c(1,2,3,4) # ROC time points
tx
predM_lp[1:10]
survROC_lp <- create_survivalROC_data(tx, dataM_test, predM_lp)

survROC_lp2 <- survROC_lp %>% mutate(sens = TP, spec = 1- FP, 
               num = sens*prev_test, 
               den = sens*prev_test + (1-spec)*(1-prev_test),
               PPV = num/den,
               Youden = TP- FP,
               ER01 = sqrt(FP**2 +(1-TP)**2),  # The Closest to (0,1) Criteria (ER) Pepe (2003), Perkins (2006)
               Liu  =sens*spec,                # Concordance Probability Method (CZ), Liu 2012
               IU  = abs(sens -auc) + abs(spec -auc) # Index of Union (IU): Unal (2017)
         ) %>% select(-c(num, den))
 
survROC_lp2 %>% print (n=50)



## ----step1-res------------------------------------------------------
step_sel
survROC_list[[step_sel]]  <- survROC_lp2
survROC_lp2 <- NULL


## ----surv-roc-plot--------------------------------------------------
## Plot Time-dependent ROC every 2.5 years
plot_timedep_ROC(survROC_lp) 


## ----step15-roc, echo = FALSE---------------------------------------
step_sel <- "step15"


## ----surv-roc15-----------------------------------------------------
dataM_test <- test_data
dataM_test$predM_lp <- NULL


lmbda_opt <-  lmbda_sel[names(lmbda_sel) == step_sel]
predM_lp <- as.vector(predict(pfit_aopt, newx = x_test, type = "link", s = lmbda_opt))
# predM_lp <- - predM_lp # !!!
predM_lp[1:10]
summary(predM_lp)

                
## Augment `dataM_test` with linear predictor
 dataM_test$predM_lp <- predM_lp

# Evaluate every 2.5 years
 
tx <- 2.5* c(1,2,3,4) # ROC time points
tx
predM_lp[1:10]
survROC_lp <- create_survivalROC_data(tx, dataM_test, predM_lp)

survROC_lp2 <- survROC_lp %>% mutate(sens = TP, spec = 1- FP, 
               num = sens*prev_test, 
               den = sens*prev_test + (1-spec)*(1-prev_test),
               PPV = num/den,
               Youden = TP- FP,
               ER01 = sqrt(FP**2 +(1-TP)**2),  # The Closest to (0,1) Criteria (ER) Pepe (2003), Perkins (2006)
               Liu  =sens*spec,                # Concordance Probability Method (CZ), Liu 2012
               IU  = abs(sens -auc) + abs(spec -auc) # Index of Union (IU): Unal (2017)
               
         ) %>% select(-c(num, den))
 
survROC_lp2 %>% print (n=50)



## ----step15-res-----------------------------------------------------
step_sel
survROC_list[[step_sel]]  <- survROC_lp2
survROC_lp2 <- NULL


## ----surv-roc15-plot------------------------------------------------
## Plot Time-dependent ROC every 2.5 years
plot_timedep_ROC(survROC_lp) 


## ----step18-roc, echo = FALSE---------------------------------------
step_sel <- "step18"


## ----surv-roc18-----------------------------------------------------
dataM_test <- test_data
dataM_test$predM_lp <- NULL


lmbda_opt <-  lmbda_sel[names(lmbda_sel) == step_sel]
predM_lp <- as.vector(predict(pfit_aopt, newx = x_test, type = "link", s = lmbda_opt))
# predM_lp <- - predM_lp # !!!
predM_lp[1:10]
summary(predM_lp)

                
## Augment `dataM_test` with linear predictor
 dataM_test$predM_lp <- predM_lp

# Evaluate every 2.5 years
 
tx <- 2.5* c(1,2,3,4) # ROC time points
tx
predM_lp[1:10]
survROC_lp <- create_survivalROC_data(tx, dataM_test, predM_lp)

survROC_lp2 <- survROC_lp %>% mutate(sens = TP, spec = 1- FP, 
               num = sens*prev_test, 
               den = sens*prev_test + (1-spec)*(1-prev_test),
               PPV = num/den,
               Youden = TP- FP,
               ER01 = sqrt(FP**2 +(1-TP)**2),  # The Closest to (0,1) Criteria (ER) Pepe (2003), Perkins (2006)
               Liu  =sens*spec,                # Concordance Probability Method (CZ), Liu 2012
               IU  = abs(sens -auc) + abs(spec -auc) # Index of Union (IU): Unal (2017)
         ) %>% select(-c(num, den))
 
survROC_lp2 %>% print (n=50)



## ----step18-res-----------------------------------------------------
step_sel
survROC_list[[step_sel]]  <- survROC_lp2
survROC_lp2 <- NULL


## ----surv-roc18-plot------------------------------------------------
## Plot Time-dependent ROC every 2.5 years
plot_timedep_ROC(survROC_lp) 


## ----step24-roc, echo = FALSE---------------------------------------
step_sel <- "step24"


## ----surv-roc24-----------------------------------------------------
dataM_test <- test_data
dataM_test$predM_lp <- NULL

lmbda_opt <-  lmbda_sel[names(lmbda_sel) == step_sel]
predM_lp <- as.vector(predict(pfit_aopt, newx = x_test, type = "link", s = lmbda_opt))
# predM_lp <- - predM_lp # !!!
predM_lp[1:10]
summary(predM_lp)

                
## Augment `dataM_test` with linear predictor
 dataM_test$predM_lp <- predM_lp

# Evaluate every 2.5 years
 
tx <- 2.5* c(1,2,3,4) # ROC time points
tx
predM_lp[1:10]
survROC_lp <- create_survivalROC_data(tx, dataM_test, predM_lp)

survROC_lp2 <- survROC_lp %>% mutate(sens = TP, spec = 1- FP, 
               num = sens*prev_test, 
               den = sens*prev_test + (1-spec)*(1-prev_test),
               PPV = num/den,
               Youden = TP- FP,
               ER01 = sqrt(FP**2 +(1-TP)**2),  # The Closest to (0,1) Criteria (ER) Pepe (2003), Perkins (2006)
               Liu  =sens*spec,                # Concordance Probability Method (CZ), Liu 2012
               IU  = abs(sens -auc) + abs(spec -auc) # Index of Union (IU): Unal (2017)

         ) %>% select(-c(num, den))
 
survROC_lp2 %>% print (n=50)



## ----step24-res-----------------------------------------------------
step_sel
survROC_list[[step_sel]]  <- survROC_lp2
survROC_lp2 <- NULL


## ----surv-roc24-plot------------------------------------------------
## Plot Time-dependent ROC every 2.5 years
plot_timedep_ROC(survROC_lp) 


## ----step31-roc, echo = FALSE---------------------------------------
step_sel <- "step31"


## ----surv-roc31-----------------------------------------------------
dataM_test <- test_data
dataM_test$predM_lp <- NULL
lmbda_opt <-  lmbda_sel[names(lmbda_sel) == step_sel]
predM_lp <- as.vector(predict(pfit_aopt, newx = x_test, type = "link", s = lmbda_opt))
# predM_lp <- - predM_lp # !!!
predM_lp[1:10]
summary(predM_lp)

                
## Augment `dataM_test` with linear predictor
 dataM_test$predM_lp <- predM_lp

# Evaluate every 2.5 years
 
tx <- 2.5* c(1,2,3,4) # ROC time points
tx
predM_lp[1:10]
survROC_lp <- create_survivalROC_data(tx, dataM_test, predM_lp)
survROC_lp2 <- survROC_lp %>% mutate(sens = TP, spec = 1- FP, 
               num = sens*prev_test, 
               den = sens*prev_test + (1-spec)*(1-prev_test),
               PPV = num/den,
               Youden = TP- FP,
               ER01 = sqrt(FP**2 +(1-TP)**2),  # The Closest to (0,1) Criteria (ER) Pepe (2003), Perkins (2006)
               Liu  =sens*spec,                # Concordance Probability Method (CZ), Liu 2012
               IU  = abs(sens -auc) + abs(spec -auc) # Index of Union (IU): Unal (2017)

         ) %>% select(-c(num, den))
 
survROC_lp2 %>% print (n=50)



## ----step31-res-----------------------------------------------------
step_sel
survROC_list[[step_sel]]  <- survROC_lp2
survROC_lp2 <- NULL


## ----surv-roc-plot31------------------------------------------------
## Plot Time-dependent ROC every 2.5 years
plot_timedep_ROC(survROC_lp) 


## ----save-xlsx------------------------------------------------------
library(writexl)
str(survROC_list)
xlsxf <- "survROC_list.xlsx"
xlsxp <- paste0("./validation_save/", xlsxf)
write_xlsx(survROC_list, xlsxp)



## ----c_index_diff-function------------------------------------------
library(boot)

c_index_diff <- function(data, indices) {
  # Create the bootstrapped dataset
  data_boot <- data[indices,]
  x_boot <- x_test[indices,] # x
  
  # Get the survival object for the bootstrapped data
  surv_obj_boot <- with(data_boot, Surv(time, status))
  
  # Generate predictions for each model on this bootstrapped dataset
  pred_lambda_min_boot <- predict(cvfit, newx = x_boot, s = 'lambda.min', type = 'response')
  pred_lambda_1se_boot <- predict(cvfit, newx = x_boot, s = 'lambda.1se', type = 'response')
  
  # Calculate C-index for each set of predictions
  #- c_index_min_boot <- concordance(Surv(time, status) ~ pred_lambda_min_boot, data = data_boot)$concordance
  #- c_index_1se_boot <- concordance(Surv(time, status) ~ pred_lambda_1se_boot, data = data_boot)$concordance

  c_index_min_boot <- concordance(Surv(time, status) ~ pred_lambda_min_boot, data = data_boot)$concordance
  c_index_1se_boot <- concordance(Surv(time, status) ~ pred_lambda_1se_boot, data = data_boot)$concordance
  
  # Return the difference
  return(c_index_min_boot - c_index_1se_boot)
}


## ----c-index-diff-bootstrap-----------------------------------------
# Bootstrapping with the 'boot' function
set.seed(123) # For reproducibility
boot_results <- boot(test_data, statistic = c_index_diff, R = 1000) # Use more iterations (R) for more stable estimates


boot_est <- mean(boot_results$t)
boot_est 

# Calculate the confidence interval from the bootstrap results
boot_ci_perc <- boot.ci(boot_results, type = "perc") # Percentile CI; could also consider 'bca' for bias-corrected accelerated CI
boot_ci_bca <- boot.ci(boot_results, type = "bca") # Percentile CI; could also consider 'bca' for bias-corrected accelerated CI

print(boot_ci_perc)
print(boot_ci_bca)


## ----minus-2ll------------------------------------------------------
# Extract the index of lambda.min and lambda.1se from the cv.fit object
lambda_min_index <- which(cvfit$lambda == cvfit$lambda.min)
lambda_1se_index <- which(cvfit$lambda == cvfit$lambda.1se)

# Extract the penalized log partial likelihood for lambda.min and lambda.1se
log_likelihood_min <- cvfit$cvm[lambda_min_index]
log_likelihood_1se <- cvfit$cvm[lambda_1se_index]

# Compute -2 Log-likelihood
minus_2LL_min <- -2 * log_likelihood_min
minus_2LL_1se <- -2 * log_likelihood_1se

# Print the values
cat("-2 Log Likelihood for lambda.min model:", minus_2LL_min, "\n")
cat("-2 Log Likelihood for lambda.1se model:", minus_2LL_1se, "\n")


## ----AICc-compare---------------------------------------------------
# Extract coefficients and calculate degrees of freedom for lambda.min and lambda.1se

coef_min <- coef(cvfit, s = "lambda.min", exact = TRUE)[,1]
coef_1se <- coef(cvfit, s = "lambda.1se", exact = TRUE)[,1]

coef_min
coef_1se

linear_predictors_min <- as.matrix(x_test) %*% coef_min ## [-1]  # Drop the intercept term if present
linear_predictors_1se <- as.matrix(x_test) %*% coef_1se

test_data$lp_min <- linear_predictors_min
test_data$lp_1se <- linear_predictors_1se

cox_model_min <- coxph(formula = ySurv_test ~ 1 + offset(lp_min), data = test_data, x = TRUE, y = TRUE)
cox_model_1se <- coxph(formula = ySurv_test ~ 1 + offset(lp_1se), data = test_data, x = TRUE, y = TRUE)

                       
df_min <- sum(coef_min != 0)
df_1se <- sum(coef_1se != 0)

# Calculate the log partial likelihood (without the penalty term) for lambda.min and lambda.1se
loglik_min <- logLik(cox_model_min)
loglik_1se <- logLik(cox_model_1se)

# Sample size
n <- nrow(x_test)

# Compute AIC and AICc using effective degrees of freedom
aic_min <- -2 * loglik_min + 2 * df_min
aic_1se <- -2 * loglik_1se + 2 * df_1se

aicc_min <- aic_min + (2 * df_min * (df_min + 1)) / (n - df_min - 1)
aicc_1se <- aic_1se + (2 * df_1se * (df_1se + 1)) / (n - df_1se - 1)

# Print the AICc values
cat("AICc for lambda.min model:", aicc_min, "\n")
cat("AICc for lambda.1se model:", aicc_1se, "\n")


## ----exit0, include=FALSE-------------------------------------------
knitr::knit_exit()

