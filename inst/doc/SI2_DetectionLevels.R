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

tox_list <- create_toxEval(full_path)

ACClong <- get_ACC(tox_list$chem_info$CAS)
ACClong <- remove_flags(ACClong)

cleaned_ep <- clean_endPoint_info(endPointInfo)
filtered_ep <- filter_groups(cleaned_ep)

# Substitute max LDL or MDL for actual values:

tox_list$chem_data <- tox_list$chem_data %>%
  left_join(select(tox_list$chem_info,
                   CAS,
                   MDL = `Maximum method detection level`,
                   LDL = `Maximum laboratory reporting level`),
            by="CAS") %>%
  rowwise() %>%
  mutate(Value = max(MDL, LDL, na.rm = TRUE)) %>%
  select(SiteID, `Sample Date`, CAS, Value) %>%
  distinct()

chemicalSummary <- get_chemical_summary(tox_list, ACClong, filtered_ep)

#Trim some names:
levels(chemicalSummary$Class)[levels(chemicalSummary$Class) == "Antimicrobial Disinfectants"] <- "Antimicrobial"
levels(chemicalSummary$Class)[levels(chemicalSummary$Class) == "Detergent Metabolites"] <- "Detergent"
levels(chemicalSummary$Class)[levels(chemicalSummary$Class) == "Flavors and Fragrances"] <- "Flavor/Fragrance"


plot_DL <- plot_tox_boxplots(chemicalSummary, "Chemical")

plot_DL


