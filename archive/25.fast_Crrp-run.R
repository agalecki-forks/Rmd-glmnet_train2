 options(width = 70)
 
# Clears Global environment
 rm(list=ls())

 bnm <-"25.fast_Crrp"   # !!! Basename for Rmd file
 mod_lbl = "M3"          # !!! do not change it
 survSplit_cut <- 10     # Time horizon 10 years
 
 output_postfix ="M3"     # !!!
 
 # user name
  username <- Sys.getenv("USER")
  if (username == "") username <- Sys.getenv("USERNAME")
 nmRmd <- paste0(bnm, ".Rmd")
 nmR <- paste0("./purl/", bnm,".Rprog")
 knitr::purl(nmRmd, output = nmR)
 output_file <- paste0(bnm,"_", output_postfix)
 rmarkdown::render(nmRmd, "all", output_file = output_file) 
