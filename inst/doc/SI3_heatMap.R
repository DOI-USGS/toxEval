## ----setup, include=FALSE---------------------------------
library(knitr)
library(rmarkdown)
options(continue=" ")
options(width=60)
knitr::opts_chunk$set(echo = TRUE,
                      warning = FALSE,
                      message = FALSE,
                      fig.height = 7,
                      fig.width = 7)

## ---------------------------------------------------------
library(readxl)
library(toxEval)
library(dplyr)
library(tidyr)
library(ggplot2)

path_to_tox <-  system.file("extdata", package="toxEval")
file_name <- "OWC_data_fromSup.xlsx"
full_path <- file.path(path_to_tox, file_name)

chem_data <- read_excel(full_path, sheet = "Data")
chem_info <- read_excel(full_path, sheet = "Chemicals") 
chem_site <- read_excel(full_path, sheet = "Sites")

#Remove DEET:
chem_data <- chem_data[chem_data$CAS != "134-62-3",] 
chem_info <- chem_info[chem_info$CAS != "134-62-3",] 
#Trim names for graph:
chem_info$Class[chem_info$Class == "Human Drug, Non Prescription"] <- "Human Drug"
chem_info$Class[chem_info$Class == "Antimicrobial Disinfectant"] <- "Antimicrobial"
chem_info$Class[chem_info$Class == "Detergent Metabolites"] <- "Detergent"

ACClong <- get_ACC(chem_info$CAS)
ACClong <- remove_flags(ACClong)

cleaned_ep <- clean_endPoint_info(endPointInfo)
filtered_ep <- filter_groups(cleaned_ep)

chemicalSummary <- get_chemical_summary(ACClong,
                                        filtered_ep,
                                        chem_data, 
                                        chem_site, 
                                        chem_info)

# Order the Great Lakes:
chem_site$site_grouping <- factor(chem_site$site_grouping,
               levels=c("Lake Superior",
               "Lake Michigan",
               "Lake Huron",
               "Lake Erie",
               "Lake Ontario"))

# Order sites:
 sitesOrdered <- c("StLouis","Nemadji","WhiteWI","Bad","Montreal",
                   "PresqueIsle","Ontonagon","Sturgeon","Tahquamenon","Burns",
                   "IndianaHC","StJoseph","PawPaw","Kalamazoo","GrandMI",
                   "Milwaukee","Muskegon","WhiteMI","PereMarquette","Manitowoc",
                   "Manistee","Fox","Oconto","Peshtigo","Menominee",
                   "Indian","Cheboygan","Ford","Escanaba","Manistique",
                   "ThunderBay","AuSable","Rifle","Saginaw","BlackMI",
                   "Clinton","Rouge","HuronMI","Raisin","Maumee",
                   "Portage","Sandusky","HuronOH","Vermilion","BlackOH",
                   "Rocky","Cuyahoga","GrandOH","Cattaraugus","Tonawanda",
                   "Genesee","Oswego","BlackNY","Oswegatchie","Grass",
                   "Raquette","StRegis")

 chem_site$`Short Name` <- factor(chem_site$`Short Name`,
               levels = sitesOrdered)

   
 plot_tox_heatmap(chemicalSummary,
              chem_site,
              category = "Chemical")


