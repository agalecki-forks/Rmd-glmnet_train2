 options(width = 70)
 
 # Clears Global environment
 rm(list=ls())
 
  source("./Rfun/create_survivalROC_data.R")
  source("./Rfun/survivalROC_helper.R")

names2_aux <- function(basenm, row_no){ 
   # row_no: 1 through 6
   nmRmd <- paste0(basenm, ".Rmd")
   nmR <- paste0("./purl/", basenm, "_", row_no, ".Rprog")
   nm_out <- paste0(bnm,"_", row_no)
   c(nmRmd = nmRmd, nmR = nmR, nm_out = nm_out)  
}

#---<<<--- 
bnm <- "55.Cox_tidy_valid"      #  Basename of Rmd file (do not change it)
parms <- list(row_no = 1)       #  Row_no (select one out of): 1,2,3,4,5,6 

#---<<<--- Do not make changes below
nms <- names2_aux(bnm, parms$row_no)
knitr::purl(nms["nmRmd"], output = nms["nmR"])
rmarkdown::render(nms["nmRmd"], "all", output_file = nms["nm_out"], params = parms)
 
 
 
