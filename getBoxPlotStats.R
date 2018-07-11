library(readxl)
library(toxEval)

path_to_tox <-  "D:/LADData/RCode/toxEval"
file_name <- "WolfCreek.xlsx"

full_path <- file.path(path_to_tox, file_name)

chem_data <- read_excel(full_path, sheet = "Data")
chem_info <- read_excel(full_path, sheet = "Chemicals") 
chem_site <- read_excel(full_path, sheet = "Sites")

ACClong <- get_ACC(chem_info$CAS)
ACClong <- remove_flags(ACClong)

cleaned_ep <- clean_endPoint_info(endPointInfo)
filtered_ep <- filter_groups(cleaned_ep, groupCol = "intended_target_family")

chemicalSummary <- get_chemical_summary(ACClong,
                                        filtered_ep,
                                        chem_data, 
                                        chem_site, 
                                        chem_info)
bioBoxPlot <- plot_tox_boxplots(chemicalSummary, 
                              category = "Biological",
                              plot_ND = FALSE)

bioBoxPlot

gD <- graphData(chemicalSummary)
boxplotStats <- boxplot(meanEAR ~ category, data = gD)
bpStats <- boxplotStats[["stats"]]

bpStats <- as.data.frame(bpStats)
names(bpStats) <- boxplotStats[["names"]]

rownames(bpStats) <- c("lower whisker",
                       "lower hinge",
                       "median",
                       "upper hinge",
                       "extreme upper whisker")

transpose_bp <- data.frame(t(bpStats))
