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

path_to_tox <-  system.file("extdata", package="toxEval")
file_name <- "OWC_data_fromSup.xlsx"
full_path <- file.path(path_to_tox, file_name)

chem_data <- read_excel(full_path, sheet = "Data")
chem_info <- read_excel(full_path, sheet = "Chemicals") 
chem_site <- read_excel(full_path, sheet = "Sites")

#Remove DEET:
chem_data <- chem_data[chem_data$CAS != "134-62-3",] 
chem_info <- chem_info[chem_info$CAS != "134-62-3",] 

#Trim names and order for graph:
chem_info$Class[chem_info$Class == "Human Drug, Non Prescription"] <- "Human Drug"
chem_info$Class[chem_info$Class == "Antimicrobial Disinfectant"] <- "Antimicrobial"
chem_info$Class[chem_info$Class == "Detergent Metabolites"] <- "Detergent"

chem_info$Class <- factor(chem_info$Class) 
chem_info$Class <- factor(chem_info$Class, 
                          levels =  c("Detergent","Antioxidant","Herbicide",
                          "Plasticizer","Insecticide","Human Drug",
                          "Antimicrobial","Fire Retardant","PAH",
                          "Flavor/Fragrance","Solvent","Dye/pigment","Fuel",
                          "Other"))

ACClong <- get_ACC(chem_info$CAS)
ACClong <- remove_flags(ACClong)
cleaned_ep <- clean_endPoint_info(endPointInfo)
filtered_ep <- filter_groups(cleaned_ep)


## ----warning=FALSE, message=FALSE-------------------------

library(dplyr)
library(tidyr)
library(DT)

chem_info_SI <- chem_info %>%
  select(`OWC Class`=Class, 
         `Compound Name` = `Chemical Name`,
         `CAS Registry Number` = CAS,
         `Molecular Weight` = MlWt,
         EEF_avg_in.vitro,
         EEF_max_in.vitro_or_in.vivo,
         AqT_EPA_acute,
         AqT_EPA_chronic,
         AqT_other_acute) %>%
  mutate(Units = "ug/l")
         
# write.csv(chem_info_SI, file="chem_info_SI.csv", row.names = FALSE, na = "")

chem_info_SI <- chem_info_SI %>%
  mutate(EEF_avg_in.vitro = as.numeric(EEF_avg_in.vitro),
         EEF_max_in.vitro_or_in.vivo = as.numeric(EEF_max_in.vitro_or_in.vivo),
         AqT_EPA_acute = as.numeric(AqT_EPA_acute),
         AqT_EPA_chronic = as.numeric(AqT_EPA_chronic),
         AqT_other_acute = as.numeric(AqT_other_acute))
  
datatable(chem_info_SI,  rownames = FALSE) %>%
  formatRound(columns=c('EEF_max_in.vitro_or_in.vivo',
                     'EEF_avg_in.vitro',
                     'AqT_EPA_acute',
                     'AqT_EPA_chronic',
                     'AqT_other_acute'), digits=3)


## ---------------------------------------------------------

intended_target <- select(endPointInfo, intended_target_family, intended_target_family_sub, endPoint = assay_component_endpoint_name, source = assay_source_long_name) %>%
  right_join(select(filtered_ep, endPoint), by = "endPoint") %>%
  arrange(intended_target_family, intended_target_family_sub) 

intended_target$intended_target_family_sub["Zebrafish" == intended_target$intended_target_family] <- "Zebrafish"

intended_target <- intended_target %>%
  rename(`Intended Target Family`=intended_target_family,
         `Intended Target Family Sub-Family` = intended_target_family_sub) %>%
  data.frame()

datatable(intended_target,  rownames = FALSE)
# write.csv(intended_target, file="intended_target.csv", row.names = FALSE, na = "") 

## ---------------------------------------------------------

ACC <- ACC

ACC_OWC <- ACC %>%
filter


## ---------------------------------------------------------

file_name <- "AOP_crosswalk.csv"
full_path <- file.path(path_to_tox, file_name)

AOP_crosswalk <- read.csv(full_path, stringsAsFactors = FALSE)

AOP_crosswalk <- select(AOP_crosswalk, 
                        gene_symbol=Target.Gene.Symbol, 
                        AOP=AOP.name)
  
AOP <- data.frame(gene_symbol = unique(AOP_crosswalk$gene_symbol),
                  AOP = "",
                  stringsAsFactors = FALSE)
for(gene in AOP$gene_symbol){
  AOP$AOP[AOP$gene_symbol %in% gene] <- paste(AOP_crosswalk$AOP[AOP_crosswalk$gene_symbol %in% gene],collapse = ", ")
}

chemicalSummary <- get_chemical_summary(ACClong,
                                        filtered_ep,
                                        chem_data, 
                                        chem_site, 
                                        chem_info)

chemicalSummary <- chemicalSummary %>%
  left_join(select(endPointInfo, 
                   endPoint=assay_component_endpoint_name,
                   subFamily=intended_target_family_sub,
                   gene_symbol=intended_target_gene_symbol), by="endPoint") %>%
  left_join(select(chem_info, casrn=CAS, `Chemical Name`), by="casrn")
  
tableData <- chemicalSummary %>%
  rename(Chemical=`Chemical Name`) %>%
  group_by(site, endPoint, Family=Bio_category, subFamily, gene_symbol, date, Chemical) %>%
  summarize(sumEAR = sum(EAR)) %>%
  group_by(site, endPoint, Family, subFamily, gene_symbol, Chemical) %>%
  summarize(maxEAR = max(sumEAR)) %>%
  group_by(endPoint, Family, subFamily, gene_symbol, Chemical) %>%
  summarize(nSites = sum(maxEAR > 10^-3)) %>%
  data.frame() %>%
  filter(nSites > 0) %>%
  spread(Chemical, nSites) %>%
  arrange(Family, subFamily, gene_symbol) %>%
  select(endPoint, Family, subFamily, gene_symbol, everything())

chem_info <- chem_info %>%
  left_join(distinct(select(chemicalSummary, casrn, `Chemical Name`)), by=c("CAS"="casrn", "Chemical Name")) %>%
  arrange(Class, `Chemical Name`)


tableData2 <- select(tableData, -endPoint, -Family, -subFamily, -gene_symbol)
tableData$nChems <- apply(tableData2, MARGIN = 1, function(x) sum(x>0, na.rm = TRUE))

orderedCols <- chem_info$`Chemical Name`[chem_info$`Chemical Name` %in% names(tableData)]

tableData <- tableData %>%
  left_join(AOP,by="gene_symbol")

tableData <- tableData[,c("Family", "subFamily","gene_symbol",
                          "endPoint","AOP", "nChems",orderedCols)] 


datatable(tableData,  rownames = FALSE)
# write.csv(tableData, file="wholeEnchilada.csv", row.names = FALSE, na = "")

