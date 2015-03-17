## ----setup, include=FALSE---------------------------------
library(xtable)
options(continue=" ")
options(width=60)
library(knitr)


## ---------------------------------------------------------

library("webchem")

inchk <- cts_convert(query = '3380-34-5', from = 'CAS', to = 'inchikey')
info <- cts_compinfo(inchikey = inchk)
info$molweight


## ---------------------------------------------------------
library(toxEval)
filePath <- system.file("extdata", package="toxEval")
passiveData <- load(file.path(filePath, "passiveData.RData"))

head(passiveData)


