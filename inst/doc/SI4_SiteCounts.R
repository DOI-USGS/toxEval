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
exclusion <- read_excel(full_path, sheet = "Exclude")

#Trim names for graph:
chem_info$Class[chem_info$Class == "Antimicrobial Disinfectants"] <- "Antimicrobial"
chem_info$Class[chem_info$Class == "Detergent Metabolites"] <- "Detergent"
chem_info$Class[chem_info$Class == "Flavors and Fragrances"] <- "Flavor/Fragrance"

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

tableData <- chemicalSummary %>%
  group_by(site, date, chnm) %>% 
  summarize(sumEAR = sum(EAR)) %>%
  group_by(site, chnm) %>%
  summarize(meanEAR = max(sumEAR)) %>%
  group_by(chnm) %>%
  summarize(nSites = sum(meanEAR > 10^-3)) %>%
  data.frame() %>%
  arrange(desc(nSites)) %>%
  filter(nSites > 1)

tableData$chnm <- factor(tableData$chnm, levels = tableData$chnm)

chemPlot <- ggplot(tableData)+
  geom_bar(aes(x=chnm, y=nSites),stat = "identity",fill = "steelblue") +
  theme_bw() +
  xlab("") +
  ylab("Number of Sites\n with EARmax > 0.001") +
  theme(axis.text.x = element_text( angle = 90,vjust=0.5,hjust = 1)) 

chemPlot


