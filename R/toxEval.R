.onAttach <- function(libname, pkgname) {
  packageStartupMessage("This information is preliminary or provisional and
is subject to revision. It is being provided to meet
the need for timely best science. The information
has not received final approval by the U.S. Geological
Survey (USGS) and is provided on the condition that
neither the USGS nor the U.S. Government shall be held
liable for any damages resulting from the authorized
or unauthorized use of the information.
                        
Funding for toxEval expires summer 2018, 
after which bug fixes & new features will be minimal")
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

#' Constants included with toxEval. Units are log micro-Molarity (log uM).
#' 
#' AC50gain downloaded on October 2015 from ToxCast dashboard. AC50gain values 
#' are either the reported modl_ga (winning model) or 10% of modl_ga if the AC50gain
#' value is lower than the lowest measured concentration. Also, 
#' 
#'\itemize{
#'  \item{ACC}{ACC endpoints}
#'}
#'
#'@aliases ACC
#'@name Constants
#'@docType data
#'@export ACC
#'@keywords datasets
#'@examples
#'ACCColumnNames <- names(ACC)
NULL

#' Endpoint information from ToxCast
#' 
#' Downloaded on October 2015 from ToxCast dashboard
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

#' tox_chemicals 
#' 
#' Downloaded on October 2015 from ToxCast dashboard
#'
#'@aliases tox_chemicals
#'@name tox_chemicals
#'@docType data
#'@keywords datasets
#'@export tox_chemicals
#'@examples
#'head(tox_chemicals)
NULL