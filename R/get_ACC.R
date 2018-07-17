#' Get the ACC values for a selection of chemicals
#' 
#' The \code{get_ACC} function retrieves the activity concentration at cutoff 
#' (ACC) values for specified chemicals. 
#' 
#' The data used in toxEval were combined from files in the 
#' "INVITRODB_V2_LEVEL5" directory that were included in the October 2015 
#' release of the ToxCast database. The function \code{get_ACC} will 
#' convert the ACC values in the ToxCast database from units of (log \eqn{\mu}M) 
#' to units of \eqn{\mu}g/L, and reformat the data as input to toxEval.
#' 
#' @param CAS Vector of CAS. 
#' @return data frame with columns CAS, chnm, flags, endPoint, ACC, MlWt, and ACC_value
#' @export
#' @importFrom tidyr gather
#' @importFrom dplyr select filter right_join mutate
#' @examples
#' CAS <- c("121-00-6","136-85-6","80-05-7","84-65-1","5436-43-1","126-73-8")
#' ACC <- get_ACC(CAS)
get_ACC <- function(CAS){
  
  # Getting rid of NSE warnings:
  Structure_MolWt <- Substance_CASRN <- casn <- chnm <- flags <- MlWt <- ACC_value <- casrn <- endPoint <- ".dplyr"
  
  chem_list <- select(tox_chemicals,
                    casrn = Substance_CASRN,
                    MlWt = Structure_MolWt) %>%
    filter(casrn %in% CAS) 
  
  ACC <- ACC %>%
    filter(CAS %in% CAS) %>%
    right_join(chem_list,
               by= c("CAS"="casrn")) %>%
    data.frame() %>%
    mutate(ACC_value = 10^ACC) %>%
    mutate(ACC_value = ACC_value * MlWt) %>%
    filter(!is.na(ACC_value)) 

  if(any(is.na(ACC$MlWt))){
    warning("Some chemicals are missing molecular weights")
  }
  
  return(ACC)

}