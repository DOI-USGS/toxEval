## ----end_point_info, eval=FALSE-----------------------------------------------
#  library(toxEval)
#  end_point_info <- end_point_info

## ----eval=FALSE---------------------------------------------------------------
#  filtered_ep <- filter_groups(end_point_info,
#                groupCol = "intended_target_family",
#                assays = c("ATG","NVS", "OT", "TOX21",
#                           "CEETOX", "APR", "CLD", "TANGUAY",
#                           "NHEERL_PADILLA","NCCT_SIMMONS", "ACEA"),
#                remove_groups = c("Background Measurement",
#                                  "Undefined"))

## ---- eval=FALSE--------------------------------------------------------------
#  install.packages("toxEval")

## ----update, echo=FALSE-------------------------------------------------------
knitr::include_graphics("update.png")

## ----eval=FALSE---------------------------------------------------------------
#  update.packages()

## ---- eval=TRUE---------------------------------------------------------------
citation(package = "toxEval")

