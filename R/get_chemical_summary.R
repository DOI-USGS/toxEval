#' get_chemical_summary
#' 
#' Get ACC values for vector of CAS's
#' @param ACClong data frame with columns: casn, chnm, MlWt, endPoint, ACC_value
#' @param filtered_ep data frame with colums: endPoints, groupCol
#' @param chem.data data frame with (at least) columns: CAS, SiteID, Value
#' @param chem.site data frame with (at least) columns: SiteID, Short Name
#' @param chem.info data frame with (at least) columns: CAS, class
#' @param exclusion data frame with columns "CAS" and "endPoint"
#' @export
#' @importFrom tidyr gather
#' @importFrom dplyr full_join filter mutate select left_join right_join anti_join
#' @examples
#' library(readxl)
#' path_to_tox <-  system.file("extdata", package="toxEval")
#' file_name <- "OWC_data_fromSup.xlsx"
#' full_path <- file.path(path_to_tox, file_name)
#' 
#' chem_data <- read_excel(full_path, sheet = "Data")
#' chem_info <- read_excel(full_path, sheet = "Chemicals") 
#' chem_site <- read_excel(full_path, sheet = "Sites")
#' 
#' ACClong <- get_ACC(chem_info$CAS)
#' ACClong <- remove_flags(ACClong)
#' 
#' cleaned_ep <- clean_endPoint_info(endPointInfo)
#' filtered_ep <- filter_groups(cleaned_ep)
#' 
#' chemicalSummary <- get_chemical_summary(ACClong,
#'                                         filtered_ep,
#'                                        chem_data, 
#'                                         chem_site, 
#'                                         chem_info)
#'                                         
#' exclusion <- read_excel(full_path, sheet = "Exclude")
#' chemicalSummary1 <- get_chemical_summary(ACClong,
#'                                         filtered_ep,
#'                                        chem_data, 
#'                                         chem_site, 
#'                                         chem_info,
#'                                         exclusion)
get_chemical_summary <- function(ACClong, filtered_ep,
                                 chem.data, chem.site, chem.info,exclusion=NULL){

  # Getting rid of NSE warnings:
  casn <- chnm <- MlWt <- endPoint <- ACC_value <- Value <- `Sample Date` <- SiteID <- ".dplyr"
  EAR <- `Short Name` <- CAS <- Class <- site <- casrn <- groupCol <- ".dplyr"
  
  if(class(chem.data$Value) == "character"){
    chem.data$Value <- as.numeric(chem.data$Value)
  }
  
  chemicalSummary <- full_join(select(ACClong, 
                                      casn, chnm, MlWt, endPoint, ACC_value), 
                               chem.data[,c("CAS", "SiteID", "Value", "Sample Date")], by=c("casn"="CAS")) %>%
    filter(!is.na(ACC_value)) %>%
    filter(!is.na(Value)) %>%
    mutate(EAR = Value/ACC_value) %>%
    rename(site = SiteID,
           date = `Sample Date`,
           casrn = casn) %>%
    select(casrn, chnm, endPoint, site, date, EAR) %>%
    filter(endPoint %in% filtered_ep$endPoint) %>%
    left_join(chem.site[,c("SiteID", "Short Name")],
              by=c("site"="SiteID")) %>%
    left_join(chem.info[, c("CAS", "Class")], by=c("casrn"="CAS")) %>%
    left_join(select(filtered_ep, endPoint, groupCol), by="endPoint") %>%
    rename(Bio_category = groupCol,
           shortName = `Short Name`)
  
  if(!is.null(exclusion)){
    chemicalSummary <- exclude_points(chemicalSummary, exclusion)
  }
  
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


#' Exclude endPoint/Chem combos
#' 
#' Using a dataframe "exclusion", filter out all chemical/endpoint combos
#' 
#' @param chemicalSummary data frame from \code{graph_chem_data}
#' @param exclusion data frame with columns "CAS" and "endPoint"
#' 
#' @export
#' @importFrom dplyr filter
#' @importFrom dplyr anti_join
#' @examples 
#' library(readxl)
#' path_to_tox <-  system.file("extdata", package="toxEval")
#' file_name <- "OWC_data_fromSup.xlsx"
#' full_path <- file.path(path_to_tox, file_name)
#' 
#' chem_data <- read_excel(full_path, sheet = "Data")
#' chem_info <- read_excel(full_path, sheet = "Chemicals") 
#' chem_site <- read_excel(full_path, sheet = "Sites")
#' exclusion <- read_excel(full_path, sheet = "Exclude")
#' ACClong <- get_ACC(chem_info$CAS)
#' ACClong <- remove_flags(ACClong)
#' 
#' cleaned_ep <- clean_endPoint_info(endPointInfo)
#' filtered_ep <- filter_groups(cleaned_ep)
#' 
#' chemicalSummary <- get_chemical_summary(ACClong,
#'                                         filtered_ep,
#'                                        chem_data, 
#'                                         chem_site, 
#'                                         chem_info,
#'                                         exclusion)
#' chemicalSummary <- exclude_points(chemicalSummary, exclusion)
exclude_points <- function(chemicalSummary, exclusion){
  
  CAS <- endPoint <- casrn <- ".dplyr"
  
  exclude_chem <- exclusion$CAS[is.na(exclusion$endPoint)]
  exclude_ep <- exclusion$endPoint[is.na(exclusion$CAS)]
  
  exclude_combo <- exclusion %>%
    filter(!is.na(CAS),
           !is.na(endPoint))
  
  chem_filtered <- chemicalSummary %>%
    filter(!(casrn %in% exclude_chem),
           !(endPoint %in% exclude_ep)) 
  
  if(nrow(exclude_combo) > 0){
    chem_filtered <- chem_filtered %>%
      anti_join(exclude_combo, by=c("casrn"="CAS",
                                    "endPoint"))
  }

  
  return(chem_filtered)
}
