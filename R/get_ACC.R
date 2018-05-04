#' Get the ACC values for a selection of chemicals
#' 
#' The get_ACC function will retrieve the activity concentration at 
#' cutoff (ACC) values for the specified chemicals. The data was originally 
#' downloaded for toxEval on October 2015 from ToxCast. The data were 
#' combined from files in the "INVITRODB_V2_LEVEL5" directory. At the time 
#' of the toxEval package release, this information was found here in the 
#' "ToxCast & Tox21 Data Spreadsheet" data set. The function get_ACC will 
#' convert the ACC values in the ToxCast database from units of (log \eqn{\mu}M) 
#' to units of \eqn{\mu}g/L, and reformat the data to a format that can be used 
#' in toxEval.
#' 
#' @param CAS vector of CAS. 
#' @return data frame with columns CAS, chnm, flags, endPoint, ACC, MlWt, and ACC_value
#' @export
#' @importFrom tidyr gather
#' @importFrom dplyr select filter right_join mutate
#' @examples
#' CAS <- c("121-00-6","136-85-6","80-05-7","84-65-1","5436-43-1","126-73-8")
#' ACClong <- get_ACC(CAS)
get_ACC <- function(CAS){
  
  # Getting rid of NSE warnings:
  Structure_MolWt <- Substance_CASRN <- casn <- chnm <- flags <- MlWt <- ACC_value <- casrn <- endPoint <- ".dplyr"
  
  chem_list <- select(tox_chemicals,
                    casrn = Substance_CASRN,
                    MlWt = Structure_MolWt) %>%
    filter(casrn %in% CAS) 
  
  ACClong <- ACC %>%
    filter(casn %in% CAS) %>%
    gather(endPoint, ACC, -casn, -chnm, -flags) %>%
    right_join(chem_list,
               by= c("casn"="casrn")) %>%
    data.frame() %>%
    mutate(ACC_value = 10^ACC) %>%
    mutate(ACC_value = ACC_value * MlWt) %>%
    filter(!is.na(ACC_value)) %>%
    rename(CAS = casn)

  if(any(is.na(ACClong$MlWt))){
    warning("Some chemicals are missing molecular weights")
  }
  
  return(ACClong)

}