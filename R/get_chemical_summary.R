#' get_chemical_summary
#' 
#' Get ACC values for vector of CAS's
#' @param ACClong data frame with columns: casn, chnm, MlWt, endPoint, ACC_value
#' @param filtered_ep data frame with colums: endPoints, groupCol
#' @param chem.data data frame with (at least) columns: CAS, SiteID, Value
#' @param chem.site data frame with (at least) columns: SiteID, Short Name
#' @param chem.info data frame with (at least) columns: CAS, class
#' @export
#' @importFrom tidyr gather
#' @importFrom dplyr full_join filter mutate select left_join right_join
#' @examples
#' CAS <- c("121-00-6","136-85-6","80-05-7","84-65-1","5436-43-1","126-73-8")
#' ACClong <- get_ACC(CAS)
#' ACClong <- remove_flags(ACClong)
#' cleaned_ep <- clean_endPoint_info(endPointInfo)
#' filtered_ep <- filter_groups(cleaned_ep)
#' 
get_chemical_summary <- function(ACClong, filtered_ep,
                                 chem.data, chem.site, chem.info){

  # Getting rid of NSE warnings:
  casn <- chnm <- MlWt <- endPoint <- ACC_value <- Value <- `Sample Date` <- SiteID <- ".dplyr"
  EAR <- `Short Name` <- CAS <- Class <- site <- casrn <- ".dplyr"
  
  chemicalSummary <- full_join(select(ACClong, 
                                      casn, chnm, MlWt, endPoint, ACC_value), 
                               chem.data, by=c("casn"="CAS")) %>%
    filter(!is.na(ACC_value)) %>%
    filter(!is.na(Value)) %>%
    mutate(EAR = Value/ACC_value) %>%
    rename(site = SiteID,
           date = `Sample Date`,
           casrn = casn) %>%
    select(casrn, chnm, endPoint, site, date, EAR) %>%
    filter(endPoint %in% filtered_ep$endPoint) %>%
    left_join(select(chem.site, SiteID, shortName = `Short Name`),
              by=c("site"="SiteID")) %>%
    left_join(select(chem.info, CAS, Class), by=c("casrn"="CAS")) %>%
    left_join(select(filtered_ep, endPoint, groupCol), by="endPoint") %>%
    rename(Bio_category = groupCol)
  
  chemicalSummary <- 

  return(chemicalSummary)
}

#' remove_flags
#' 
#' Remove endpoints with specific flags from data
#' 
#' @param ACClong data frame with columns: casn, chnm, MlWt, endPoint, ACC_value
#' @param flagsShort vector of flags to TAKE OUT. Possible values are 
#' "Borderline", "OnlyHighest", "OneAbove","Noisy", "HitCall", "GainAC50", "Biochemical"
#' @export
#' @examples 
#' CAS <- c("121-00-6","136-85-6","80-05-7","84-65-1","5436-43-1","126-73-8")
#' ACClong <- get_ACC(CAS)
#' ACClong <- remove_flags(ACClong)
remove_flags <- function(ACClong, flagsShort = c("Borderline",
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
  
  flags <- ".dplyr"
  
  flag_hits <- select(ACClong, flags) %>%
    mutate(Borderline = grepl("Borderline active", flags),
           Noisy = grepl("Noisy data", flags),
           OneAbove = grepl("Only one conc above baseline", flags),
           OnlyHighest = grepl("Only highest conc above baseline", flags),
           Biochemical = grepl("Biochemical assay with", flags),
           GainAC50 = grepl("Gain AC50", flags),
           HitCall = grepl("potentially confounded by overfitting", flags)) %>%
    select(-flags)
  
  ACClong <- ACClong[rowSums(flag_hits[flagsShort]) == 0,]

  return(ACClong)
  
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