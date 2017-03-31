#' filter_groups
#' 
#' Clean up the endPointInfo table from toxCast. Filtering and cleaning based on ES&T (cite Dan/Brett's paper)
#' 
#' @param ep data frame Endpoint information from ToxCast
#' @param groupCol character name of column to use as a group catetory
#' @param assays vector of assays to use. Possible values are "ATG","NVS","OT","TOX21","CEETOX", "APR", "BSK",
#' "CLD","TANGUAY","NHEERL_PADILLA","NCCT_SIMMONS","ACEA" 
#' @param remove_groups vector of groups to remove
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
                          remove_groups = "Background Measurement"){
  
  match.arg(assays, c("ATG","NVS","OT","TOX21","CEETOX", "APR", "BSK",
                      "CLD","TANGUAY","NHEERL_PADILLA",
                      "NCCT_SIMMONS","ACEA"), several.ok = TRUE)
  
  ep <- ep[,c("assay_component_endpoint_name",groupCol,"assay_source_name")] %>%
    rename(endPoint = assay_component_endpoint_name,
           assaysFull = assay_source_name)
  names(ep)[names(ep) == groupCol] <- "groupCol"
  
  ep <- ep[(ep$assaysFull %in% assays),]
  ep <- ep[!is.na(ep$groupCol),]
  ep <- ep[!(ep$groupCol) %in% remove_groups,]
  
  return(ep)
}



#' filter_flags
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
#' 
filter_flags <- function(ep, flagsShort = c("Borderline",
                                            "OnlyHighest",
                                            "GainAC50",
                                            "Biochemical")){

  match.arg(flagsShort, 
            c("Borderline",
              "OnlyHighest",
              "OneAbove",
              "Noisy",
              "HitCall",
              "GainAC50",
              "Biochemical"),
            several.ok = TRUE)
  # So, with the defaults, we are taking out:
  # c("Borderline active",
  #   "Only highest conc above baseline, active",
  #   "Gain AC50 < lowest conc & loss AC50 < mean conc", 
  #   "Biochemical assay with < 50% efficacy")
  # We are leaving in with the defaults:
  # c("Hit-call potentially confounded by overfitting",
  #   "Only one conc above baseline, active",
  #   "Noisy data")
}