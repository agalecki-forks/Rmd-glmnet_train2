 options(width = 70)
 
# Clears Global environment
 rm(list=ls())

 bnm <-"20.create_npx"   # !!! Basename for Rmd file
 survSplit_cut <- 10     # Time horizon 10 years ( do not change it)
 
  # user name
  username <- Sys.getenv("USER")
  if (username == "") username <- Sys.getenv("USERNAME")
  
  output_postfix = substr(username,1,3) #"M3"     # !!!
  mod_lbl = "M3"          # !!! do not change it


 nmRmd <- paste0(bnm, ".Rmd")
 nmR <- paste0("./purl/", bnm,".Rprog")
 knitr::purl(nmRmd, output = nmR)
 output_file <- paste0(bnm,"_", output_postfix)
 rmarkdown::render(nmRmd, "all", output_file = output_file) 
