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
library(toxEval)
library(grid)

path_to_tox <-  system.file("extdata", package="toxEval")
file_name <- "OWC_data_fromSup.xlsx"
full_path <- file.path(path_to_tox, file_name)

tox_list <- create_toxEval(full_path)

ACClong <- get_ACC(tox_list$chem_info$CAS)
ACClong <- remove_flags(ACClong)

cleaned_ep <- clean_endPoint_info(endPointInfo)
filtered_ep <- filter_groups(cleaned_ep)

chemicalSummary <- get_chemical_summary(tox_list, ACClong, filtered_ep)

#Trim some names:
levels(chemicalSummary$Class)[levels(chemicalSummary$Class) == "Antimicrobial Disinfectants"] <- "Antimicrobial"
levels(chemicalSummary$Class)[levels(chemicalSummary$Class) == "Detergent Metabolites"] <- "Detergent"
levels(chemicalSummary$Class)[levels(chemicalSummary$Class) == "Flavors and Fragrances"] <- "Flavor/Fragrance"

bioPlot <- plot_tox_boxplots(chemicalSummary, 
                               category = "Biological", 
                               manual_remove = c("Transferase","Undefined"))

bioPlot
grid.text("# Samples", 
         x = unit(.22, "npc"), 
         y = unit(.995, "npc"), gp=gpar(fontsize=7))

## ---------------------------------------------------------
classPlot <- plot_tox_boxplots(chemicalSummary, 
                               category = "Chemical Class")
classPlot
grid.text("# Samples", 
         x = unit(.22, "npc"), 
         y = unit(.995, "npc"), gp=gpar(fontsize=7))

## ---------------------------------------------------------

chemPlot <- plot_tox_boxplots(chemicalSummary, 
                               category = "Chemical")
chemPlot
grid.text("# Samples", 
         x = unit(.35, "npc"), 
         y = unit(.995, "npc"), gp=gpar(fontsize=7))

