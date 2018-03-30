## ------------------------------------------------------------------------
library(toxEval)
path_to_tox <-  system.file("extdata", package="toxEval")
file_name <- "OWC_data_fromSup.xlsx"
full_path <- file.path(path_to_tox, file_name)

tox_list <- create_toxEval(full_path)


## ----chemicalSummary-----------------------------------------------------


