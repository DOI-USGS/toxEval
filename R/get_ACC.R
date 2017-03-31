#' get_ACC
#' 
#' Get ACC values for vector of CAS's
#' 
#' @param CAS vector of CAS
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
    filter(!is.na(ACC_value)) 

  
  return(ACClong)

}