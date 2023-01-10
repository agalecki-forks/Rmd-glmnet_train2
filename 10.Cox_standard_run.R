 options(width = 70)
 
 # Clears Global environment
 rm(list=ls())
 
 names_aux <- function(basenm, mod_lbl){ 
   # mod_lbl M0, M1 or M2
   nmRmd <- paste0(basenm, ".Rmd")
   nmR <- paste0("./purl/", basenm, "_", mod_lbl, ".Rprog")
   nm_out <- paste0(bnm,"_", mod_lbl)
   c(nmRmd = nmRmd, nmR = nmR, nm_out = nm_out)  
}

#---<<<--- 
bnm <- "10.Cox_standard"            #  Basename of Rmd file
parms <- list(mod_lbl = "M2")        #  Model lbl: M0, M1, M2 

#---<<<--- Do not make changes below
nms <- names_aux(bnm, parms$mod_lbl)
knitr::purl(nms["nmRmd"], output = nms["nmR"])
rmarkdown::render(nms["nmRmd"], "all", output_file = nms["nm_out"], params = parms)
 
 
 
