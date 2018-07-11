library(toxEval)
library(readxl)
library(dplyr)
library(ggplot2)

path_to_tox <-  "D:/LADData/RCode/toxEval"
file_name <- "neonic_bench.xlsx"

full_path <- file.path(path_to_tox, file_name)

chem_data <- read_excel(full_path, sheet = "Data")
chem_info <- read_excel(full_path, sheet = "Chemicals") 
chem_site <- read_excel(full_path, sheet = "Sites")
exclusion <- read_excel(full_path, sheet = "Exclude")

benchmarks <- read_excel(full_path, sheet = "Benchmarks") %>%
  rename(chnm = Compound,
         ACC_value = Value) %>%
  filter(!is.na(CAS))

filtered_ep <- select(benchmarks, endPoint) %>%
  distinct() %>%
  mutate(groupCol = "Aquatic Benchmark")

chemicalSummary <- get_chemical_summary(benchmarks,
                                        filtered_ep,
                                        chem_data, 
                                        chem_site, 
                                        chem_info,
                                        exclusion)

chem_plot <- plot_tox_boxplots(chemicalSummary, category = "Chemical")

ggsave(chem_plot, file = "aquatic_benchmarks.pdf", height = 11, width=9)

ep_plot <- plot_tox_endpoints(chemicalSummary)
ggsave(ep_plot, file = "aquatic_benchmarks_by_endpoint.pdf", height = 5, width=7)

table_tox_endpoint(chemicalSummary, category = "Chemical")
plot_tox_stacks(chemicalSummary, chem_site, category = "Chemical", include_legend = FALSE)

