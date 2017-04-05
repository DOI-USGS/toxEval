## ----setup, include=FALSE---------------------------------
library(knitr)
library(rmarkdown)
options(continue=" ")
options(width=60)
knitr::opts_chunk$set(echo = TRUE)

## ----warning=FALSE, message=FALSE, fig.width=6, fig.height=6----
library(readxl)
library(toxEval)

path_to_tox <-  system.file("extdata", package="toxEval")
file_name <- "OWC_data_fromSup.xlsx"

full_path <- file.path(path_to_tox, file_name)

chem_data <- read_excel(full_path, sheet = "Data")
chem_info <- read_excel(full_path, sheet = "Chemicals") 
chem_site <- read_excel(full_path, sheet = "Sites")

#Remove DEET:
chem_data <- chem_data[chem_data$CAS != "134-62-3",] 
chem_info <- chem_info[chem_info$CAS != "134-62-3",] 

ACClong <- get_ACC(chem_info$CAS)
ACClong <- remove_flags(ACClong)

cleaned_ep <- clean_endPoint_info(endPointInfo)
filtered_ep <- filter_groups(cleaned_ep)

chemicalSummary <- get_chemical_summary(ACClong,
                                        filtered_ep,
                                        chem_data, 
                                        chem_site, 
                                        chem_info)

#Trim some names:
chemicalSummary$Class[chemicalSummary$Class == "Human Drug, Non Prescription"] <- "Human Drug"
chemicalSummary$Class[chemicalSummary$Class == "Antimicrobial Disinfectant"] <- "Antimicrobial"
chemicalSummary$Class[chemicalSummary$Class == "Detergent Metabolites"] <- "Detergent"

bioPlot <- plot_group_boxplots(chemicalSummary, 
                               filtered_ep, 
                               category = "Biological", 
                               manual_remove = c("Transferase","Undefined"))
bioPlot

## ----warning=FALSE, fig.width=6, fig.height=6-------------
classPlot <- plot_group_boxplots(chemicalSummary,
                                 filtered_ep,
                                 category = "Chemical Class")
classPlot


## ----warning=FALSE, fig.width=6, fig.height=6-------------
chemPlot <- plot_chemical_boxplots(chemicalSummary,
                                 filtered_ep)
chemPlot


