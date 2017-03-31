#' get_chemical_summary
#' 
#' Get ACC values for vector of CAS's
#' 
#' @param ACClong data frame with columns: casn, chnm, MlWt, endPoint, ACC_value
#' @param filtered_ep vector of end points
#' @param chem.data data frame with (at least) columns: CAS, SiteID, Value
#' @param chem.site data frame with (at least) columns: SiteID, Short Name
#' @param chem.info data frame with (at least) columns: CAS, class
#' @param flagsShort vector of flags to TAKE OUT. Possible values are 
#' "Borderline", "OnlyHighest", "OneAbove","Noisy", "HitCall", "GainAC50", "Biochemical"
#' @export
#' @importFrom tidyr gather
#' @importFrom dplyr full_join filter mutate select left_join right_join
#' @examples
#' CAS <- c("121-00-6","136-85-6","80-05-7","84-65-1","5436-43-1","126-73-8")
#' ACClong <- get_ACC(CAS)
#' cleaned_ep <- clean_endPoint_info(endPointInfo)
#' filtered_ep <- filter_groups(cleaned_ep)
#' 
get_chemical_summary <- function(ACClong, filtered_ep,
                                 chem.data, chem.site, chem.info){
  
  match.arg(flagsShort,
            c("Borderline", "OnlyHighest", "OneAbove","Noisy", "HitCall", "GainAC50", "Biochemical"),
            several.ok = TRUE)
  
  casn <- chnm <- MlWt <- endPoint <- ACC_value <- Value <- `Sample Date` <- SiteID <- ".dplyr"
  
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
    filter(endPoint %in% filtered_ep) %>%
    left_join(select(chem.site, SiteID, shortName = `Short Name`),
              by=c("site"="SiteID"))

  chemicalSummary <- chemicalSummary %>%
    left_join(select(chem.info, CAS, Class), by=c("casrn"="CAS"))

  return(chemicalSummary)
}


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
  
  for(i in flagsShort){
    take.out.flags <- flagDF[!flagDF[[i]],c("casn","endPoint")]

    chemicalSummary <- right_join(chemicalSummary, take.out.flags,
                                  by=c("casrn"="casn", "endPoint"="endPoint")) %>%
      filter(!is.na(chnm))
  }
  
}