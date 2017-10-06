## ----setup, include=FALSE---------------------------------
library(knitr)
library(rmarkdown)
options(continue=" ")
options(width=60)
knitr::opts_chunk$set(echo = TRUE,
                      warning = FALSE,
                      message = FALSE,
                      fig.height = 7,
                      fig.width = 10)

## ---------------------------------------------------------
library(readxl)
library(toxEval)
library(dplyr)
library(ggplot2)
library(grid)

path_to_tox <-  system.file("extdata", package="toxEval")
file_name <- "OWC_data_fromSup.xlsx"

full_path <- file.path(path_to_tox, file_name)

chem_data <- read_excel(full_path, sheet = "Data")
chem_info <- read_excel(full_path, sheet = "Chemicals") 
chem_site <- read_excel(full_path, sheet = "Sites")
exclusion <- read_excel(full_path, sheet = "Exclude")

ACClong <- get_ACC(chem_info$CAS)
ACClong <- remove_flags(ACClong)

cleaned_ep <- clean_endPoint_info(endPointInfo)
filtered_ep <- filter_groups(cleaned_ep)

chemicalSummary <- get_chemical_summary(ACClong,
                                        filtered_ep,
                                        chem_data, 
                                        chem_site, 
                                        chem_info,
                                        exclusion)

#Trim some names:
chemicalSummary$Class[chemicalSummary$Class == "Human Drug, Non Prescription"] <- "Human Drug"
chemicalSummary$Class[chemicalSummary$Class == "Antimicrobial Disinfectant"] <- "Antimicrobial"
chemicalSummary$Class[chemicalSummary$Class == "Detergent Metabolites"] <- "Detergent"

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

chem_site <- chem_site %>%
  mutate(`Short Name` = factor(`Short Name`, levels = sitesOrdered)) %>%
  mutate(site_grouping = factor(site_grouping, levels = c("Lake Superior",
                                                          "Lake Michigan",
                                                          "Lake Huron",
                                                          "Lake Erie",
                                                          "Lake Ontario")))

for(class in unique(chemicalSummary$Class)){
  sub_class <- filter(chemicalSummary, Class %in% class)
  
  upperPlot <- plot_tox_stacks(sub_class, 
                               chem_site, 
                               category = "Chemical")
  upperPlot <- upperPlot +
    ggtitle(class)
    
  print(upperPlot)
  grid.text("# Samples:", 
          x = unit(.03, "npc"), 
          y = unit(.205, "npc"), gp=gpar(fontsize=7))
}




