## ----end_point_info, eval=FALSE------------------------------------------
#  library(toxEval)
#  end_point_info <- end_point_info

## ----eval=FALSE----------------------------------------------------------
#  filtered_ep <- filter_groups(end_point_info,
#                groupCol = "intended_target_family",
#                assays = c("ATG","NVS", "OT", "TOX21",
#                           "CEETOX", "APR", "CLD", "TANGUAY",
#                           "NHEERL_PADILLA","NCCT_SIMMONS", "ACEA"),
#                remove_groups = c("Background Measurement",
#                                  "Undefined"))

## ---- eval=FALSE---------------------------------------------------------
#  rprofile_path = file.path(Sys.getenv("HOME"), ".Rprofile")
#  write('\noptions(repos=c(getOption(\'repos\'),
#      CRAN=\'https://cloud.r-project.org\',
#      USGS=\'https://owi.usgs.gov/R\'))\n',
#        rprofile_path,
#        append =  TRUE)
#  
#  cat('Your Rprofile has been updated to include GRAN.
#      Please restart R for changes to take effect.')

## ---- eval=FALSE---------------------------------------------------------
#  install.packages("toxEval")

## ----eval=FALSE----------------------------------------------------------
#  update.packages()

## ---- eval=FALSE---------------------------------------------------------
#  citation(package = "toxEval")

