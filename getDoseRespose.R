library(toxEval)
library(tcpl)
library(RMySQL)
library(data.table)
library(readxl)
library(dplyr)

tcplConf(drvr = "MySQL", 
         user = "toxEval",
         pass = "mysql",
         host = "130.11.179.26",
         db = "toxcastdb")

ToxCastInfoPath <- "D:/LADData/RCode/toxEval_Archive/INVITRODB_V2_SUMMARY"
ChemSummary <- "Chemical_Summary_151020.csv"
m4id <- "m4id_Matrix_151020.csv"

dfChemSum <- read.csv(file.path(ToxCastInfoPath, ChemSummary),stringsAsFactors = FALSE)
dfm4id <- read.csv(file.path(ToxCastInfoPath, m4id),stringsAsFactors = FALSE)
names(dfm4id)[1] <- "ChemCode"
# 
# #Get chemical code for Atrizine
# chCode <- dfChemSum[["code"]][which(dfChemSum$casn %in% "1912-24-9")]
# 
# #Determine which row represents caffiene
# chemRow <- which(dfm4id[,"ChemCode"] %in% chCode)
# 
# m4id <- dfm4id[chemRow,"TOX21_p53_BLA_p1_ratio"]
# 
# tcplPlotM4ID(m4id = m4id, lvl = 5)

path_to_file <- 'D:/LADData/RCode/toxEval_Archive/Scripts for prepping data/new_tw2016.xlsx' 
tox_list <- create_toxEval(path_to_file)
ACClong <- get_ACC(tox_list$chem_info$CAS)
ACClong <- remove_flags(ACClong = ACClong,
                        flagsShort = c('Borderline','OnlyHighest','GainAC50','Biochemical'))

cleaned_ep <- clean_endPoint_info(endPointInfo)
filtered_ep <- filter_groups(cleaned_ep, 
                             groupCol = 'intended_target_family',
                             assays = c('ATG','NVS','OT','TOX21','CEETOX','APR','CLD','TANGUAY','NHEERL_PADILLA','NCCT_SIMMONS','ACEA'),
                             remove_groups = c('Background Measurement','Undefined'))

chemicalSummary <- get_chemical_summary(tox_list, ACClong, filtered_ep)

chem_names <- select(chemicalSummary, CAS, chnm) %>%
  distinct()

owc_cas <- create_toxEval("D:/LADData/RCode/toxEval_Archive/Scripts for Paper/OWC_data_fromSup.xlsx")
owc_cas <- owc_cas$chem_info %>%
  pull(CAS)

pest_cas <- create_toxEval("D:/LADData/RCode/GLRI_CEC_2016/data/pesticides.xlsx")
pest_cas <- pest_cas$chem_info %>%
  pull(CAS)

chem_names <- select(chemicalSummary, CAS, chnm) %>%
  distinct() %>%
  filter(!(CAS %in% c(pest_cas, owc_cas)))

order_chms <- chemicalSummary %>%
  filter(CAS %in% chem_names$CAS) %>%
  group_by(chnm, CAS) %>%
  summarise(sumEAR = sum(EAR, na.rm = TRUE)) %>%
  arrange(desc(sumEAR))

for(i in 1:nrow(order_chms)){
  
  chnm <- as.character(order_chms$chnm[i])
  chem <- order_chms$CAS[i]
  
  sub_df <- filter(chemicalSummary, CAS == chem) %>%
    select(endPoint, EAR) %>%
    distinct() %>%
    group_by(endPoint) %>%
    summarise(EAR = max(EAR)) %>%
    arrange(desc(EAR))
  
  chCode <- dfChemSum[["code"]][which(dfChemSum$casn %in% chem)]
  chemRow <- which(dfm4id[,"ChemCode"] %in% chCode)
  
  if(length(sub_df[["endPoint"]]) > 0){
  
    pdf(paste0("D:/LADData/RCode/toxEval_Archive/Scripts for prepping data/pBradley_dose/",chnm,".pdf"), width = 11, height = 8)
    for(ep_to_plot in sub_df[["endPoint"]]){
      m4id <- dfm4id[chemRow,ep_to_plot]
      
      tryCatch({
        tcplPlotM4ID(m4id = m4id, lvl = 5)
      }, error = function(e){
        message("no data",ep_to_plot)
      })
    }
    dev.off()
  }
  
}
####################################
suspect <- read.csv("D:/LADData/RCode/GLRI_CEC_2016/suspect_endpoints.csv", stringsAsFactors = FALSE, na.strings = "")

names(suspect) <- c("CAS","endPoint")

exclude_chem <- suspect$CAS[is.na(suspect$endPoint)]
exclude_ep <- suspect$endPoint[is.na(suspect$CAS)]

exclude_combo <- suspect %>%
  filter(!is.na(CAS),
         !is.na(endPoint))

pathToPests <- "D:/LADData/RCode/GLRI_CEC_2016/data"
file_name <- "pesticides.xlsx"

full_path <- file.path(pathToPests, file_name)

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

pdf("D:/LADData/RCode/GLRI_CEC_2016/dose_suspect.pdf", width = 11, height = 8)

for(i in 1:nrow(exclude_combo)){
  
  chCode <- dfChemSum[["code"]][which(dfChemSum$casn %in% exclude_combo$CAS[i])]
  chemRow <- which(dfm4id[,"ChemCode"] %in% chCode)
  ep_to_plot <- exclude_combo$endPoint[i]
  m4id <- dfm4id[chemRow,ep_to_plot]
    
  tryCatch({
    tcplPlotM4ID(m4id = m4id, lvl = 5)
  }, error = function(e){
    message("no data",ep_to_plot)
  })
  
}

dev.off()

pdf("D:/LADData/RCode/GLRI_CEC_2016/dose_suspect_full_chems.pdf", width = 11, height = 8)

for(i in exclude_chem){
  
  chCode <- dfChemSum[["code"]][which(dfChemSum$casn %in% i)]
  chemRow <- which(dfm4id[,"ChemCode"] %in% chCode)
  
  sub_df <- filter(chemicalSummary, CAS == i) %>%
    select(endPoint, EAR) %>%
    distinct() %>%
    group_by(endPoint) %>%
    summarise(EAR = max(EAR)) %>%
    arrange(desc(EAR))
  
  
  for(ep_to_plot in sub_df[["endPoint"]]){
    m4id <- dfm4id[chemRow,ep_to_plot]
    
    tryCatch({
      tcplPlotM4ID(m4id = m4id, lvl = 5)
    }, error = function(e){
      message("no data",ep_to_plot)
    })
  }
  
}

dev.off()

###############################################

suspect <- read.csv("D:/LADData/RCode/GLRI_CEC_2016/10_load_data/exclusions.csv", stringsAsFactors = FALSE, na.strings = "")

names(suspect) <- c("CAS","endPoint")

exclude_chem <- suspect$CAS[is.na(suspect$endPoint)]
exclude_ep <- suspect$endPoint[is.na(suspect$CAS)]

exclude_combo <- suspect %>%
  filter(!is.na(CAS),
         !is.na(endPoint))


pdf("D:/LADData/RCode/GLRI_CEC_2016/dose_excluded.pdf", width = 11, height = 8)
for(i in 1:nrow(exclude_combo)){
  
  chCode <- dfChemSum[["code"]][which(dfChemSum$casn %in% exclude_combo$CAS[i])]
  chemRow <- which(dfm4id[,"ChemCode"] %in% chCode)
  ep_to_plot <- exclude_combo$endPoint[i]
  m4id <- dfm4id[chemRow,ep_to_plot]
  
  tryCatch({
    tcplPlotM4ID(m4id = m4id, lvl = 5)
  }, error = function(e){
    message("no data",ep_to_plot)
  })
  
}

dev.off()

pdf("D:/LADData/RCode/GLRI_CEC_2016/dose_excluded_full_eps.pdf", width = 11, height = 8)
for(i in exclude_ep){

  m4id <- dfm4id[chemRow,ep_to_plot]

  chCode <- dfChemSum[["code"]][which(dfChemSum$casn %in% unique(chemicalSummary$CAS))]
  chemRow <- which(dfm4id[,"ChemCode"] %in% chCode)
    
  for(j in chemRow){
    m4id <- dfm4id[j,i]
    tryCatch({
      tcplPlotM4ID(m4id = m4id, lvl = 5)
    }, error = function(e){
      message("no data",ep_to_plot)
    })
  }
}

dev.off()
