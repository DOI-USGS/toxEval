library(toxEval)
library(dplyr)

#Change these next 2 lines:
path_file <-  system.file("extdata", package="toxEval")
file_name <- "OWC_data_fromSup.xlsx"

full_path <- file.path(path_file, file_name)

tox_list <- create_toxEval(full_path)

ACClong <- get_ACC(tox_list$chem_info$CAS)
ACClong <- remove_flags(ACClong)

cleaned_ep <- clean_endPoint_info(endPointInfo)
filtered_ep <- filter_groups(cleaned_ep)

chemicalSummary <- get_chemical_summary(tox_list, ACClong, filtered_ep)
endPointInfo <- endPointInfo

chemicalSummary <- chemicalSummary %>%
  left_join(select(endPointInfo, 
                   endPoint=assay_component_endpoint_name, 
                   intended_target_family,
                   intended_target_family_sub, 
                   intended_target_official_symbol), by="endPoint")
