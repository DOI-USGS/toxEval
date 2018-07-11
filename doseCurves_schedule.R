library(toxEval)
library(tcpl)
library(RMySQL)
library(data.table)
library(readxl)
library(dplyr)
library(tidyr)

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

schedule <- "2437" #pesticide
url.info <- paste0("http://nwqlqc/servlets_u/AnalyticalServicesGuide?srchStr=",schedule,"&srchType=sched&mCrit=exact&oFmt=xl")
temp.path <- tempdir()
temp.file <- file.path(temp.path,"AnalyticalServicesGuide.xls")
download.file(url.info, destfile = temp.file, mode="wb")

schedule.data <- read_excel(temp.file)
ACC <- ACC

cas_in_tox <- unique(ACC$casn)
cas_in_tox_schedule <- cas_in_tox[which(cas_in_tox %in% unique(schedule.data$`CAS Number`))]

ACC_sch <- ACC %>%
  filter(casn %in% cas_in_tox_schedule) %>%
  gather(endPoint, value, -chnm, -flags, -casn) %>%
  filter(!is.na(value)) %>%
  arrange(casn, value)

for(i in seq_along(cas_in_tox_schedule)){
  
  cas <- cas_in_tox_schedule[i]
  chnm <- schedule.data$`Parameter Name`[schedule.data$`CAS Number` == cas]
  endpoints <- filter(ACC_sch, casn == cas)[["endPoint"]]
  
  chCode <- dfChemSum[["code"]][which(dfChemSum$casn %in% cas)]
  chemRow <- which(dfm4id[,"ChemCode"] %in% chCode)
  
  if(length(endpoints) > 0){
    pdf(paste0("D:/LADData/RCode/toxEval_Archive/Does Response/Pesticide/",chnm,".pdf"), width = 11, height = 8)
    for(ep_to_plot in endpoints){
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
