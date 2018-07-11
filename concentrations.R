library(readxl)
library(toxEval)
library(dplyr)
library(tidyr)

path_to_tox <-  "D:/LADData/RCode/toxEval"
file_name <- "OWC_data_fromSup_concentration.xlsx"
full_path <- file.path(path_to_tox, file_name)

chem_data <- read_excel(full_path, sheet = "Data")
chem_info <- read_excel(full_path, sheet = "Chemicals") 
chem_site <- read_excel(full_path, sheet = "Sites")
exclusion <- read_excel(full_path, sheet = "Exclude")

#Trim names and order for graph:
chem_info$Class[chem_info$Class == "Antimicrobial Disinfectants"] <- "Antimicrobial"
chem_info$Class[chem_info$Class == "Detergent Metabolites"] <- "Detergent"
chem_info$Class[chem_info$Class == "Flavors and Fragrances"] <- "Flavor/Fragrance"

chem_info$Class <- factor(chem_info$Class) 

ACClong <- get_ACC(chem_info$CAS)
ACClong <- remove_flags(ACClong)
cleaned_ep <- clean_endPoint_info(endPointInfo)
filtered_ep <- filter_groups(cleaned_ep)

benchmarks <- read_excel(full_path, sheet = "Benchmarks") %>%
  rename(chnm = `Chemical Name`,
         ACC_value = value) %>%
  filter(!is.na(CAS))

chemicalSummary <- get_chemical_summary(benchmarks,
                                        filtered_ep,
                                        chem_data, 
                                        chem_site, 
                                        chem_info,
                                        exclusion)



filtered_ep <- select(benchmarks, endPoint) %>%
  distinct() %>%
  mutate(groupCol = "Concentrations")

chem_plot <- plot_tox_boxplots(chemicalSummary, category = "Chemical")
chem_plot
