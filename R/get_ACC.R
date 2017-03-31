#' filter_groups
#' 
#' Clean up the endPointInfo table from toxCast. Filtering and cleaning based on ES&T (cite Dan/Brett's paper)
#' 
#' @param ep data frame Endpoint information from ToxCast
#' @param flagsShort vector of flags TO REMOVE
#' @export
#' @importFrom stringi stri_trans_totitle
#' @examples 
#' endPointInfo <- endPointInfo
#' cleaned_ep <- clean_endPoint_info(endPointInfo)
#' filtered_ep <- filter_groups(cleaned_ep)
filter_groups <- function(ep, 
                          groupCol = "intended_target_family",
                          assays = c("ATG","NVS","OT","TOX21","CEETOX","APR", 
                                     "CLD","TANGUAY","NHEERL_PADILLA",
                                     "NCCT_SIMMONS","ACEA"),
                          remove_group = "Background Measurement"){
  
  match.arg(assays, c("ATG","NVS","OT","TOX21","CEETOX", "APR", "BSK",
                      "CLD","TANGUAY","NHEERL_PADILLA",
                      "NCCT_SIMMONS","ACEA"), several.ok = TRUE)
  
  ep <- ep[,c("assay_component_endpoint_name",groupCol,"assay_source_name")] %>%
    rename(endPoint = assay_component_endpoint_name,
           assaysFull = assay_source_name)
  names(ep)[names(ep) == groupCol] <- "groupCol"
  
  ep <- ep[(ep$assaysFull %in% assays),]
  ep <- ep[!is.na(ep$groupCol),]
  ep <- ep[!(ep$groupCol) %in% remove_group,]
  
  return(ep)
}



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