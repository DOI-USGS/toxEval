.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Although this software program has been used by the U.S. Geological Survey (USGS), no warranty, expressed or implied, is made by the USGS or the U.S. Government as to the accuracy and functioning of the program and related program material nor shall the fact of distribution constitute any such warranty, and no responsibility is assumed by the USGS in connection therewith.")
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
#' http://www.usgs.gov/visual-id/credit_usgs.html#copyright\cr
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
#'  \item{AC50gain}{AC50 gain endpoints}
#'  \item{AC50loss}{AC50 loss endpoints}
#'  \item{AC10}{AC10 endpoints}
#'}
#'
#'@aliases AC50gain AC50loss AC10
#'@name Constants
#'@docType data
#'@export AC50gain AC50loss AC10
#'@keywords datasets
#'@examples
#'AC50GainColumns <- names(AC50gain)
#'AC50LossColumns <- names(AC50loss)
#'AC10Columns <- names(AC10)
NULL

#' Endpoint information from ToxCast
#' 
#' Downloaded on October 2015 from ToxCast dashboard
#'
#'@name endPointInfo
#'@aliases endPointInfo flagDF
#'@docType data
#'@keywords datasets
#'@export endPointInfo flagDF
#'@examples
#'endPointInfo <- endPointInfo
#'head(endPointInfo[,1:5])
#'flagDF <- flagDF
#'head(flagDF)
NULL

#' passiveData
#' 
#'
#'@aliases passiveData
#'@name passiveData
#'@docType data
#'@keywords datasets
NULL

#' wData 
#'
#'@aliases wData
#'@name wData
#'@docType data
#'@keywords datasets
NULL

#' pCodeInfo 
#'
#'@aliases pCodeInfo
#'@name pCodeInfo
#'@docType data
#'@keywords datasets
NULL