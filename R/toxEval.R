.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    paste(strwrap('This information is preliminary or provisional and
is subject to revision. It is being provided to meet
the need for timely best science. The information
has not received final approval by the U.S. Geological
Survey (USGS) and is provided on the condition that
neither the USGS nor the U.S. Government shall be held
liable for any damages resulting from the authorized
or unauthorized use of the information.

USGS Research Package: 
https://owi.usgs.gov/R/packages.html#research'),
      collapse='\n'))
}


#' Evaluation of ToxCast data with measured concentrations
#'
#' \tabular{ll}{
#' Package: \tab toxEval\cr
#' Type: \tab Package\cr
#' License: \tab Unlimited for this package, dependencies have more restrictive licensing.\cr
#' Copyright: \tab This software is in the public domain because it contains materials
#' that originally came from the United States Geological Survey, an agency of
#' the United States Department of Interior. For more information, see the
#' official USGS copyright policy at
#' https://www.usgs.gov/visual-id/credit_usgs.html#copyright\cr
#' LazyLoad: \tab yes\cr
#' }
#'
#'  Initial code for studying ToxCast data in relation to measured concentrations.
#'
#' @name toxEval-package
#' @docType package
#' @importFrom dplyr filter
#' @importFrom dplyr rename
#' @importFrom dplyr group_by
#' @importFrom dplyr summarize
#' @importFrom dplyr summarise
#' @importFrom dplyr select
#' @importFrom dplyr arrange
#' @importFrom dplyr distinct
#' @importFrom dplyr left_join
#' @importFrom dplyr right_join
#' @importFrom dplyr desc
#' @importFrom dplyr filter_
#' @importFrom dplyr rename_
#' @importFrom dplyr mutate
#' @importFrom dplyr group_by_
#' @importFrom dplyr select_
#' @importFrom dplyr mutate_
#' @author Steven Corsi \email{srcorsi@@usgs.gov}, Laura De Cicco \email{ldecicco@@usgs.gov}
#' @keywords ToxCast
NULL

#' ACC values included with toxEval. 
#' 
#' Downloaded on October 2015 from ToxCast. The data were
#' combined from files in the "INVITRODB_V2_LEVEL5" folder. 
#' At the time of toxEval package release, this information was found:
#' \url{https://www.epa.gov/chemical-research/toxicity-forecaster-toxcasttm-data}
#' in the "ToxCast & Tox21 Data Spreadsheet" data set. 
#' 
#' The data has been provided in a "wide" format, however
#' the \code{get_ACC} function is an easy way to get the data
#' in a "long" format. AC50gain values are the reported modl_ga (winning model) and units are 
#' log micro-Molarity (log \eqn{\mu}M).
#' 
#'
#'@aliases ACC
#'@name ACC
#'@docType data
#'@export ACC
#'@keywords datasets
#'@examples
#'ACCColumnNames <- names(ACC)
NULL

# If we need to update the ACC data frame, here is
# is a function that *should* do it, assuming
# the format it the same. I might not count on that
# however....that is why it is an internal only function.
# This is staying commented out because it adds 
# extraneous notes:
# update_ACC <- function(path_to_files){
#   library(data.table)
#   library(dplyr)
#   library(tidyr)
#   # Data originally from:
#   # ftp://newftp.epa.gov/COMPTOX/ToxCast_Data_Oct_2015/README_INVITRODB_V2_LEVEL5.pdf
#   # https://www.epa.gov/sites/production/files/2015-08/documents/toxcast_assay_annotation_data_users_guide_20141021.pdf
#   # path_to_files <- "D:/LADData/RCode/toxEval_Archive/INVITRODB_V2_LEVEL5"
# 
#   files <- list.files(path = path_to_files)
#   
#   x <- fread(file.path(path_to_files, files[1]))
#   
#   filtered <- select(x, chnm, casn, aenm, logc_min, logc_max, modl_acc,
#                      modl, actp, modl_ga, flags, hitc,gsid_rep) %>%
#     filter(hitc == 1)
#   
#   for(i in files[-1]){
#     subX <- fread(file.path(path_to_files,i)) 
#     
#     subFiltered <- select(subX, chnm, casn, aenm, logc_min, logc_max, modl_acc,
#                           modl, actp, modl_ga, flags, hitc,gsid_rep) %>%
#       filter(hitc == 1)
#     
#     filtered <- bind_rows(filtered, subFiltered)
#   }
#   
#   ACCgain <- filter(filtered, hitc == 1) %>%
#     filter(gsid_rep == 1) %>%
#     select(casn, chnm, aenm, modl_acc, flags) %>%
#     spread(key = aenm, value = modl_acc)
#   
#   # Something we considered but decided not to do was:
#   
#   # ACCgain2 <- filter(filtered, hitc == 1) %>%
#   #   filter(gsid_rep == 1) %>%
#   #   select(casn, chnm, aenm, modl_acc, flags, logc_min) %>%
#   #   mutate(newFlag = modl_acc < logc_min) %>%
#   #   mutate(value = ifelse(newFlag, log10((10^modl_acc)/10), modl_acc)) 
#   
# }

#' Endpoint information from ToxCast
#' 
#' Downloaded on October 2015 from ToxCast. The file name of the
#' raw data was "Assay_Summary_151020.csv" from the zip file 
#' "Assay_Information_Oct_2015.zip". At the time
#' of toxEval package release, this information was found:
#' \url{https://www.epa.gov/chemical-research/toxicity-forecaster-toxcasttm-data}
#' in the section marked "Download Assay Information", in the 
#' ToxCast & Tox21 high-throughput assay information data set.
#'
#'
#'@name endPointInfo
#'@aliases endPointInfo
#'@docType data
#'@keywords datasets
#'@export endPointInfo
#'@examples
#'endPointInfo <- endPointInfo
#'head(endPointInfo[,1:5])
NULL

#' ToxCast Chemical Information 
#' 
#' Downloaded on October 2015 from ToxCast. The file name of the
#' raw data was "TOX21IDs_v4b_23Oct2014_QCdetails.xlsx", 
#' from the US EPA DSSTox DATA RELEASE OCTOBER 2015. At the time
#' of toxEval package release, this information was found:
#' \url{https://www.epa.gov/chemical-research/toxicity-forecaster-toxcasttm-data}
#' in the section marked "Download ToxCast Chemical Information". This 
#' was in the "ToxCast & Tox21 Chemicals Distributed Structure-Searchable Toxicity Database (DSSTox files)"
#' data set.
#'
#'@aliases tox_chemicals
#'@name tox_chemicals
#'@docType data
#'@keywords datasets
#'@export tox_chemicals
#'@examples
#'head(tox_chemicals)
NULL
