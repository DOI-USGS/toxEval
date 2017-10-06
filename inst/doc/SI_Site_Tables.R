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
library(DT)

path_to_tox <-  system.file("extdata", package="toxEval")
file_name <- "OWC_data_fromSup.xlsx"
full_path <- file.path(path_to_tox, file_name)

chem_data <- read_excel(full_path, sheet = "Data")
chem_info <- read_excel(full_path, sheet = "Chemicals") 
chem_site <- read_excel(full_path, sheet = "Sites")
exclusion <- read_excel(full_path, sheet = "Exclude")

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
                                        chem_info,
                                        exclusion)

chemicalSummary <- chemicalSummary %>%
  left_join(select(endPointInfo, 
                   endPoint=assay_component_endpoint_name,
                   subFamily=intended_target_family_sub,
                   gene_symbol=intended_target_gene_symbol), by="endPoint") %>%
  left_join(select(chem_info, casrn=CAS, `Chemical Name`), by="casrn")
  
chem_info <- chem_info %>%
  left_join(distinct(select(chemicalSummary, casrn, `Chemical Name`)), by=c("CAS"="casrn", "Chemical Name")) %>%
  arrange(Class, `Chemical Name`)



## ---------------------------------------------------------

tableData <- chemicalSummary %>%
    rename(Chemical=`Chemical Name`) %>%
    group_by(site, endPoint, Family=Bio_category, subFamily, gene_symbol, date, Chemical) %>%
    summarize(sumEAR = sum(EAR)) %>%
    group_by(site, endPoint, Family, subFamily, gene_symbol, Chemical) %>%
    summarize(maxEAR = max(sumEAR)) %>%
    filter(maxEAR > 0) %>%
    data.frame() %>%
    spread(Chemical, maxEAR) %>%
    arrange(site, Family, subFamily, gene_symbol) %>%
    select(site, Family, subFamily, gene_symbol,endPoint, everything()) %>%
    left_join(AOP,by="gene_symbol")

list_tables <- list()

for(i in 1:nrow(chem_site)){
  
  site <- chem_site$SiteID[i]
  site_name <- chem_site$`Short Name`[i]
  
  tableData_site <- tableData[tableData$site == site,]
    
  tableData_site <- Filter(function(x)!all(is.na(x)), tableData_site)
    
  list_tables[[2*i-1]] <- htmltools::tags$h3(site_name)
    
  if(nrow(tableData_site) > 0){
    tableData2 <- select(tableData_site, -endPoint, -Family, -subFamily, -gene_symbol, -AOP, -site)
    tableData_site$nChems <- apply(tableData2, MARGIN = 1, function(x) sum(x>0, na.rm = TRUE))
    
    orderedCols <- chem_info$`Chemical Name`[chem_info$`Chemical Name` %in% names(tableData_site)]
    
    tableData_site <- tableData_site[,c("Family", "subFamily","gene_symbol",
                              "endPoint","AOP", "nChems",orderedCols)] 
    
    dt_table <- datatable(tableData_site, rownames = FALSE, caption = "This????") %>%
                 formatSignif(columns=orderedCols, digits=3)  
    
    list_tables[[2*i]] <- dt_table
  } else {
    list_tables[[2*i]] <- htmltools::tags$h3("EAR never > 0")
  }
}



## ---------------------------------------------------------

htmltools::tagList(
  list_tables
)


