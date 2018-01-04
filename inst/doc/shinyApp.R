## ----setup, include=FALSE---------------------------------
library(knitr)
library(rmarkdown)
options(continue=" ")
options(width=60)
knitr::opts_chunk$set(echo = TRUE,
                      warning = FALSE,
                      message = FALSE)

## ---- eval=FALSE------------------------------------------
#  rprofile_path = file.path(Sys.getenv("HOME"), ".Rprofile")
#  write('\noptions(repos=c(getOption(\'repos\'),
#      CRAN=\'https://cloud.r-project.org\',
#      USGS=\'https://owi.usgs.gov/R\'))\n',
#        rprofile_path,
#        append =  TRUE)
#  
#  cat('Your Rprofile has been updated to include GRAN.
#      Please restart R for changes to take effect.')

## ---- eval=FALSE------------------------------------------
#  install.packages("toxEval")

## ----eval=FALSE-------------------------------------------
#  update.packages()

## ---- eval=FALSE------------------------------------------
#  citation(package = "toxEval")

