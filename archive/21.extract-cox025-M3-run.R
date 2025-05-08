 options(width = 70)
 
# Clears Global environment
 rm(list=ls())

 bnm <-"21.extract-cox025-M3"   # !!! Basename for Rmd file
 survSplit_cut <- 10     # Time horizon 10 years ( do not change it)
 
  # user name
  username <- Sys.getenv("USER")
  if (username == "") username <- Sys.getenv("USERNAME")
  
  output_postfix = if (username == "agalecki") "_test" else "" # 
  mod_lbl = "M3"          # !!! do not change it


 nmRmd <- paste0(bnm, ".Rmd")
 nmR <- paste0("./purl/", bnm,".Rprog")
 knitr::purl(nmRmd, output = nmR)
 output_file <- paste0(bnm, output_postfix)
 rmarkdown::render(nmRmd, "all", output_file = output_file) 
