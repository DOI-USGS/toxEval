
.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    paste(strwrap(paste("For more information:
https://doi-usgs.github.io/toxEval/
ToxCast database: version", dbVersion()), width = 40),
      collapse = "\n"
    )
  )
}

dbVersion <- function() {
  "4.1"
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
#' @name toxEval-package
#' @docType package
#' @author Laura De Cicco \email{ldecicco@@usgs.gov}. Steven Corsi
#' @keywords internal 
"_PACKAGE"

#' ACC values included with toxEval. See \code{vignette("Setting up toxEval package data", package = "toxEval")}
#' for more information on how the data was aggregated.
#'
#' Downloaded on September 2023 from ToxCast. See also
#' \url{https://www.frontiersin.org/articles/10.3389/ftox.2023.1275980/full}.
#'
#' @references U.S. EPA. 2023. ToxCast & Tox21 Summary Files. 
#' Retrieved from \url{https://www.epa.gov/chemical-research/toxicity-forecaster-toxcasttm-data}
#' on September 2023. 
#'
#' @source \url{https://www.epa.gov/comptox-tools/exploring-toxcast-data}
#'
#' @aliases ToxCast_ACC
#' @return data frame with columns CAS, chnm (chemical name), flags, endPoint, and ACC (value).
#' @name ToxCast_ACC
#' @docType data
#' @export ToxCast_ACC
#' @keywords datasets
#' @examples
#' head(ToxCast_ACC)
NULL

#' Endpoint information from ToxCast
#'
#' See \code{vignette("Setting up toxEval package data", package = "toxEval")}
#' for more information on how the data was aggregated.
#'
#'
#' @name end_point_info
#' @aliases end_point_info
#' @docType data
#' @keywords datasets
#' @references U.S. EPA. 2014. ToxCast Assay Annotation Data User Guide.
#'
#' @source \doi{10.23645/epacomptox.6062479.v3}
#' @export end_point_info
#' @return data frame with 72 columns. The columns and definitions
#' are discussed in the "ToxCast Assay Annotation Version 1.0 Data User Guide (PDF)" (see source).
#' The column "Relevance Category" was included for consideration of 
#' grouping/filtering endpoints based on user goals.
#' @examples
#' end_point_info <- end_point_info
#' head(end_point_info[, 1:5])
NULL


#' ToxCast Chemical Information
#'
#' See \code{vignette("Setting up toxEval package data", package = "toxEval")}
#' for more information on how the data was aggregated.
#' 
#' @aliases tox_chemicals
#' @name tox_chemicals
#' @return data frame 
#' @docType data
#' @keywords datasets
#' @export tox_chemicals
#' @examples
#' head(tox_chemicals)
NULL

#' ToxCast Chemical Information
#'
#' See \code{vignette("Setting up toxEval package data", package = "toxEval")}
#' for more information on how the data was aggregated.
#' 
#' @aliases flags
#' @name flags
#' @return data frame 
#' @docType data
#' @keywords datasets
#' @export flags
#' @examples
#' head(flags)
NULL

utils::globalVariables(c("CAS", "endPoint", "chnm", "flags", "site",
                         "Bio_category", "Class", "EAR",
                         "sumEAR", "value", "calc", "choice_calc", "nHits",
                         "Structure_MolWt", "casrn", "Substance_Name",
                         "MlWt", "ACC_value", "Substance_CASRN",
                         "Value", "Sample Date", "SiteID", "Short Name",
                         "groupCol", "Chemical", "logEAR", "meanEAR",
                         "median", "max_med", "choice_calc", "nHits",
                         "nSites", "Samples with hits", "nSamples", "hits",
                         "dec_lat", "dec_lon", "nSites", "name",
                         "nonZero", "maxEAR", "count", "site_grouping",
                         "index", "n", "x", "y", "max_med", "ymin", "label",
                         "ymax", "hit_label", "percentDet", "lab",
                         "aeid", "assay_component_endpoint_name", "casn", "hit_val"))

