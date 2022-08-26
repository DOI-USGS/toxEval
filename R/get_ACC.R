#' Get the ACC values for a selection of chemicals
#' 
#' The \code{get_ACC} function retrieves the activity concentration at cutoff 
#' (ACC) values for specified chemicals. 
#' 
#' The data used in toxEval were combined from files in the 
#' "INVITRODB_V3_LEVEL5" directory that were included in the October 2018 
#' release of the ToxCast database. The function \code{get_ACC} will 
#' convert the ACC values in the ToxCast database from units of (log \eqn{\mu}M) 
#' to units of \eqn{\mu}g/L, and reformat the data as input to toxEval.
#' 
#' @param CAS Vector of CAS.
#' @import dplyr
#' 
#' @return data frame with columns CAS, chnm, flags, endPoint, ACC, MlWt, and ACC_value
#' @export
#' @examples
#' CAS <- c("121-00-6","136-85-6","80-05-7","84-65-1","5436-43-1","126-73-8")
#' ACC <- get_ACC(CAS)
#' head(ACC)
get_ACC <- function(CAS){
  
  # Getting rid of NSE warnings:
  Structure_MolWt <- Substance_CASRN <- casn <- Substance_Name <- ".dplyr"
  chnm <- flags <- MlWt <- ACC_value <- casrn <- endPoint <- ".dplyr"
  
  chem_list <- dplyr::select(tox_chemicals,
                            casrn = Substance_CASRN,
                            MlWt = Structure_MolWt) 
  chem_list <- dplyr::filter(chem_list, casrn %in% CAS) 
  
  ACC <- ToxCast_ACC 
  ACC <- dplyr::filter(ACC, CAS %in% CAS)
  ACC <- dplyr::right_join(ACC, chem_list,
               by= c("CAS"="casrn")) 
  
  ACC <- mutate(ACC, ACC_value = 10^ACC,
                  ACC_value = ACC_value * MlWt) 
  ACC <- filter(ACC, !is.na(ACC_value))
  ACC <- left_join(ACC, select(tox_chemicals, 
                           CAS = Substance_CASRN, 
                           chnm = Substance_Name), by="CAS")

  if(any(is.na(ACC$MlWt))){
    warning("Some chemicals are missing molecular weights")
  }
  
  return(ACC)

}