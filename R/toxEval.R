
.onAttach <- function(libname, pkgname) {

  packageStartupMessage(
    paste(strwrap(paste('For more information:
https://rconnect.usgs.gov/toxEval_docs/
ToxCast database: version', dbVersion()), width = 40),
      collapse='\n'))
}

dbVersion <- function(){
  "3.5"
}

#' Analyze ToxCast data in relation to measured concentrations.
#' 
#' \code{toxEval} includes a set of functions to analyze, visualize, and 
#' organize measured concentration data as it relates to ToxCast data 
#' (default) or other user-selected chemical-biological interaction 
#' benchmark data such as water quality criteria. The intent of 
#' these analyses is to develop a better understanding of the potential 
#' biological relevance of environmental chemistry data. Results can 
#' be used to prioritize which chemicals at which sites may be of 
#' greatest concern. These methods are meant to be used as a screening 
#' technique to predict potential for biological influence from chemicals 
#' that ultimately need to be validated with direct biological assays. 

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
#'
#' @name toxEval-package
#' @docType package
#' @author Laura De Cicco \email{ldecicco@@usgs.gov}. Steven Corsi  
#' @keywords ToxCast
NULL

#' ACC values included with toxEval. 
#' 
#' Downloaded on January 2020 from ToxCast. The data were
#' combined from files in the "INVITRODB_V3_LEVEL5" folder. 
#' At the time of toxEval package release, this information was found:
#' \url{https://www.epa.gov/chemical-research/exploring-toxcast-data-downloadable-data}
#' in the "ToxCast & Tox21 Data Spreadsheet" data set. 
#' ACC values are the in the "ACC" column (winning model) and units are 
#' log micro-Molarity (log \eqn{\mu}M).
#' 
#' @references Toxicology, EPA's National Center for Computational (2018): ToxCast and Tox21 Data Spreadsheet. figshare. Dataset.
#'  https://doi.org/10.23645/epacomptox.6062503.v3
#'  
#' @source \url{https://www.epa.gov/chemical-research/exploring-toxcast-data-downloadable-data}
#'
#'@aliases ToxCast_ACC
#'@return data frame with columns CAS, chnm (chemical name), flags, endPoint, and ACC (value).
#'@name ToxCast_ACC
#'@docType data
#'@export ToxCast_ACC
#'@keywords datasets
#'@examples
#'head(ToxCast_ACC)
NULL
# 
# ToxCast_ACC <- readRDS("ACC_v3_5.rds")
# end_point_info <- readRDS("endpoint_info_v_35.rds")
# end_point_info$intended_target_family_sub[1968] <- "receptor tyrosine phosphatase"
# 
# tox_chemicals <- readRDS("chemicals3_5.rds")
# 
# tox_chemicals <- tox_chemicals |>
#   arrange(Structure_MolWt) |>
#   filter(!duplicated(Substance_CASRN))
# 
# save(ToxCast_ACC, end_point_info, tox_chemicals,
#      file = "R/sysdata.rda", compress = "xz")

#' Endpoint information from ToxCast
#' 
#' Downloaded on October 2018 from ToxCast. The file name of the
#' raw data was "Assay_Summary_190226.csv" from the zip file 
#' "INVITRODB_V3_1_SUMMARY" folder. At the time
#' of the toxEval package release, these data were found at:
#' \url{https://www.epa.gov/chemical-research/exploring-toxcast-data-downloadable-data}
#' in the section marked "Download Assay Information", in the 
#' ToxCast & Tox21 high-throughput assay information data set.
#'
#'
#'@name end_point_info
#'@aliases end_point_info
#'@docType data
#'@keywords datasets
#'@references U.S. EPA. 2014. ToxCast Assay Annotation Data User Guide. 
#'\url{https://www.epa.gov/chemical-research/toxcast-assay-annotation-data-user-guide}.
#'  
#'@source \url{https://doi.org/10.23645/epacomptox.6062479.v3}
#'@export end_point_info
#'@return data frame with 86 columns. The columns and definitions
#'are discussed in the "ToxCast Assay Annotation Version 1.0 Data User Guide (PDF)" (see source)
#'@examples
#'end_point_info <- end_point_info
#'head(end_point_info[,1:5])
NULL

#' ToxCast Chemical Information 
#' 
#' Downloaded on October 2015 from ToxCast. The file name of the
#' raw data was "TOX21IDs_v4b_23Oct2014_QCdetails.xlsx", 
#' from the US EPA DSSTox DATA RELEASE OCTOBER 2015. At the time
#' of toxEval package release, this information was found:
#' \url{https://www.epa.gov/chemical-research/exploring-toxcast-data-downloadable-data}
#' in the section marked "Download ToxCast Chemical Information". This 
#' was in the "ToxCast & Tox21 Chemicals Distributed Structure-Searchable Toxicity Database (DSSTox files)"
#' data set.
#'
#'@aliases tox_chemicals
#'@name tox_chemicals
#'@return data frame with columns:    
#'"Substance_Name","Substance_CASRN",
#'"Structure_MolWt" 
#'@docType data
#'@keywords datasets
#'@export tox_chemicals
#'@examples
#'head(tox_chemicals)
NULL
