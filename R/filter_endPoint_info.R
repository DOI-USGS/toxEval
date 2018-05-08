#' Filter endPoints based on groups and assays.
#' 
#' This function provides a mechanism to specify 3 levels of information in the 
#' supplied data frame 'endPointInfo' to be used in subsequent analysis steps. 
#' First, the user specifies the ToxCast assay annotation using the 'groupCol' 
#' argument, which is a column header in 'endPointInfo'. Second, the user 
#' specifies the families of assays to use. Finally, the user can choose to 
#' remove specific group(s) from the category. The default is to remove 
#' 'Background Measurement' and 'Undefined'. Choices for this should be 
#' reconsidered based on individual study objectives.
#' 
#' The default category ('groupCol') is 'intended_target_family'. Depending 
#' on the study, other categories may be more relevant. 
#' 
#' @param ep Data frame containing Endpoint information from ToxCast
#' @param groupCol Character name of ToxCast annotation column to use as a group catetory
#' @param assays Vector of assays to use in the data analysis. Possible values are "ATG","NVS","OT","TOX21","CEETOX", "APR", "BSK",
#' "CLD","TANGUAY","NHEERL_PADILLA","NCCT_SIMMONS","ACEA". By default, the 
#' "BSK" (BioSeek) assay is removed.
#' @param remove_groups Vector of groups within the selected 'groupCol' to remove
#' @export
#' @importFrom stringi stri_trans_totitle
#' @importFrom dplyr rename
#' @examples 
#' endPointInfo <- endPointInfo
#' cleaned_ep <- clean_endPoint_info(endPointInfo)
#' filtered_ep <- filter_groups(cleaned_ep)
filter_groups <- function(ep, 
                          groupCol = "intended_target_family",
                          assays = c("ATG","NVS","OT","TOX21","CEETOX","APR", 
                                     "CLD","TANGUAY","NHEERL_PADILLA",
                                     "NCCT_SIMMONS","ACEA"),
                          remove_groups = c("Background Measurement","Undefined")){
  
  match.arg(assays, c("ATG","NVS","OT","TOX21","CEETOX", "APR", "BSK",
                      "CLD","TANGUAY","NHEERL_PADILLA",
                      "NCCT_SIMMONS","ACEA"), several.ok = TRUE)
  
  # Getting rid of NSE warnings:
  assay_source_name <- assay_component_endpoint_name <- ".dplyr"
  
  ep <- ep[,c("assay_component_endpoint_name",groupCol,"assay_source_name")] %>%
    rename(endPoint = assay_component_endpoint_name,
           assaysFull = assay_source_name)
  names(ep)[names(ep) == groupCol] <- "groupCol"
  
  ep <- ep[(ep$assaysFull %in% assays),]
  ep <- ep[!is.na(ep$groupCol),]
  if(any(!is.na(remove_groups))){
    if(any(remove_groups != "")){
      ep <- ep[!(ep$groupCol) %in% remove_groups,]
    }
  }
  
  return(ep)
}

