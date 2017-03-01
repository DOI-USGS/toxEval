###########################
library(dplyr)
library(readr)
library(tidyr)

dataDir <- "D:/LADData/RCode/toxEval_Archive/INVITRODB_V2_LEVEL5"
setwd(dataDir)
files <- list.files()

x <- read_csv(files[1], 
              col_types = list(stkc = col_double()))

filtered <- select(x, chnm, casn, aenm, logc_min, logc_max, modl_acc,
                   modl, actp, modl_ga, flags, hitc,gsid_rep) 

for(i in files[c(-1)]){
  subX <- read_csv(i, col_types = list(stkc = col_double())) 
  
  subFiltered <- select(subX, chnm, casn, aenm, logc_min, logc_max, modl_acc,
                        modl, actp, modl_ga, flags, hitc,gsid_rep) 
  
  filtered <- bind_rows(filtered, subFiltered)
}

setwd("D:/LADData/RCode/toxEval")
rm(subFiltered, subX, x, dataDir, files, i)

source("getDataReady.R")

rm(ACC, ACClong, graphData, orderChem, orderClass, siteLimits, 
   take.out.flags, waterData, waterSamples, wData, wDataLong,
   assays, cleanUpNames, detColumns, filePath, flagsShort, i, newSiteKey,
   orderedLevels, packagePath, pathToApp, qualColumns, sitesOrdered, 
   siteToFind, valColumns, fancyNumbers, fancyNumbers2, simpleCap)


# pCodeInfo <- pCodeInfo %>%
#   arrange(class) %>%
#   left_join(distinct(select(ACC, casn)), by=c("casrn"="casn"))
# 
# pCodeInfo2 <- pCodeInfo %>%
#   arrange(class) %>%
#   left_join(distinct(select(filtered, casn)), by=c("casrn"="casn"))
# 
# x <- select(pCodeInfo, chnm, class, casrn)

allData <- select(filtered, -chnm) %>%
  right_join(select(pCodeInfo, casrn, chnm, parameter_nm, class, EEF_avg_in.vitro, EEF_max_in.vitro_or_in.vivo,
                    AqT_EPA_acute, AqT_EPA_chronic, AqT_other_acute, AqT_other_chronic), by=c("casn"="casrn")) %>%
  filter(casn %in% pCodeInfo$casrn) 

allSum <- allData %>%
  group_by(chnm, casn) %>%
  summarize(all_count = length(unique(aenm[!is.na(aenm)]))) %>%
  arrange(desc(all_count)) %>%
  data.frame()

repData <- allData %>%
  filter(gsid_rep == 1)

repSum <- repData %>%
  group_by(chnm) %>%
  summarize(rep_count = length(unique(aenm[!is.na(aenm)]))) %>%
  arrange(desc(rep_count))

# Now take out our stuff:

consideredData <- repData %>%
  filter(aenm %in% ep$endPoint)

consideredSum <- consideredData %>%
  group_by(chnm) %>%
  summarize(considered_count = length(unique(aenm[!is.na(aenm)]))) %>%
  arrange(desc(considered_count))

hitData <- consideredData %>%
  filter(hitc == 1)

hitSum <- hitData %>%
  group_by(chnm) %>%
  summarize(hitc_count = length(unique(aenm[!is.na(aenm)]))) %>%
  arrange(desc(hitc_count))

joinStuff <- left_join(allSum, consideredSum, by="chnm")
joinStuff <- left_join(joinStuff, hitSum, by="chnm")

flagsShort <- c("Borderline",  "OnlyHighest",
                "GainAC50", "Biochemical")
flagData <- hitData
flagDF <- flagDF
for(i in flagsShort){
  take.out.flags <- flagDF[!flagDF[[i]],c("casn","endPoint")]
  
  flagData <- right_join(flagData, take.out.flags, 
                                by=c("casn"="casn", 
                                     "aenm"="endPoint")) %>%
    filter(!is.na(chnm))
}

flagSum <- flagData %>%
  group_by(chnm) %>%
  summarize(flag_count = length(unique(aenm[!is.na(aenm)])),
            min_ACC = min(modl_acc, na.rm = TRUE)) %>%
  arrange(desc(flag_count))

joinStuff <- left_join(joinStuff, flagSum, by="chnm")

pCodeInfo$maxEEF <- pmax(pCodeInfo$EEF_avg_in.vitro, pCodeInfo$EEF_max_in.vitro_or_in.vivo, na.rm = TRUE)
pCodeInfo$minAqT <- pmin(pCodeInfo$AqT_EPA_acute, pCodeInfo$AqT_EPA_chronic, pCodeInfo$AqT_other_acute, pCodeInfo$AqT_other_acute, na.rm = TRUE)

joinStuff <- left_join(joinStuff, select(pCodeInfo, casn=casrn, class, maxEEF, minAqT, mlWt, parameter_units), by="casn")

joinStuff_final <- joinStuff %>%
  select(class, chnm, everything()) %>%
  arrange(class) %>%
  rename(`OWC class` = class,
         `Compound name` = chnm, 
         `CAS registry number` = casn,
         `Total` = all_count, 
         `Considered` = considered_count,
         `Active` = hitc_count,
         `Filtered` = flag_count) %>%
  data.frame() %>%
  mutate(min_ACC = 10^min_ACC) %>%
  mutate(min_ACC = min_ACC *mlWt) %>%
  mutate(min_ACC = signif(min_ACC, 4),
         mlWt = signif(mlWt, 4))

write.csv(joinStuff_final, "endPointCounts.csv", row.names = FALSE, na = "-")

joinStuff_final_2 <- joinStuff_final %>%
  rename(casrn=`CAS.registry.number`) %>%
  left_join(select(pCodeInfo, EEF_avg_in.vitro,EEF_max_in.vitro_or_in.vivo,casrn,
                   AqT_EPA_acute, AqT_EPA_chronic, AqT_other_acute, AqT_other_chronic)) %>%
  select(-Total, -Considered, -Active, -Filtered, -min_ACC, -minAqT, -minEEF)

write.csv(joinStuff_final_2, "chemInfo.csv", row.names = FALSE, na = "-")


###################################################
