library(toxEval)
library(dplyr)
library(data.table)
library(tidyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(stringi)

pCodeInfo <- pCodeInfo

pathToApp <- system.file("extdata", package="toxEval")
endPointInfo <- endPointInfo

endPointInfo <- endPointInfo[!(endPointInfo$assay_source_name == "ATG" & endPointInfo$signal_direction == "loss"),]
endPointInfo <- endPointInfo[!(endPointInfo$assay_source_name == "NVS" & endPointInfo$signal_direction == "gain"),]
endPointInfo <- endPointInfo[endPointInfo$assay_component_endpoint_name != "TOX21_p53_BLA_p3_ratio",]
endPointInfo <- endPointInfo[endPointInfo$assay_component_endpoint_name != "TOX21_p53_BLA_p2_viability",]

endPointInfo$intended_target_family[endPointInfo$assay_component_endpoint_name %in% 
                                      c("CLD_CYP1A1_24hr","CLD_CYP1A1_48hr","CLD_CYP1A1_6hr",
                                        "CLD_CYP1A2_24hr","CLD_CYP1A2_48hr","CLD_CYP1A2_6hr")] <- "dna binding"

endPointInfo$intended_target_family[endPointInfo$assay_component_endpoint_name %in% 
                                      c("CLD_CYP2B6_24hr","CLD_CYP2B6_48hr","CLD_CYP2B6_6hr",
                                        "CLD_CYP3A4_24hr","CLD_CYP3A4_48hr","CLD_CYP3A4_6hr",
                                        "CLD_SULT2A_48hr","CLD_UGT1A1_48hr","NVS_NR_bER",
                                        "NVS_NR_bPR","NVS_NR_cAR")] <- "nuclear receptor"

endPointInfo$intended_target_family[endPointInfo$assay_component_endpoint_name %in% 
                                      c("Tanguay_ZF_120hpf_ActivityScore",
                                        "Tanguay_ZF_120hpf_AXIS_up",
                                        "Tanguay_ZF_120hpf_BRAI_up",
                                        "Tanguay_ZF_120hpf_CFIN_up",
                                        "Tanguay_ZF_120hpf_EYE_up",
                                        "Tanguay_ZF_120hpf_JAW_up",
                                        "Tanguay_ZF_120hpf_MORT_up",
                                        "Tanguay_ZF_120hpf_OTIC_up",
                                        "Tanguay_ZF_120hpf_PE_up",
                                        "Tanguay_ZF_120hpf_PFIN_up",
                                        "Tanguay_ZF_120hpf_PIG_up",
                                        "Tanguay_ZF_120hpf_SNOU_up",
                                        "Tanguay_ZF_120hpf_SOMI_up",
                                        "Tanguay_ZF_120hpf_SWIM_up",
                                        "Tanguay_ZF_120hpf_TR_up",
                                        "Tanguay_ZF_120hpf_TRUN_up",
                                        "Tanguay_ZF_120hpf_YSE_up")] <- "zebrafish"

assays <- c("ATG","NVS","OT","TOX21","CEETOX", "APR", #"APR"?,"BSK"?
            "CLD","TANGUAY","NHEERL_PADILLA",
            "NCCT_SIMMONS","ACEA")

ep <- select(endPointInfo, 
             endPoint = assay_component_endpoint_name,
             groupCol = intended_target_family,
             assaysFull = assay_source_name) %>%
  filter(assaysFull %in% assays) %>%
  filter(!(groupCol %in% c("background measurement"))) %>% #,"cell morphology", "cell cycle"))) %>%
  filter(!is.na(groupCol))

cleanUpNames <- endPointInfo$intended_target_family
cleanUpNames <- stri_trans_totitle(cleanUpNames)
cleanUpNames[grep("Dna",cleanUpNames)] <- "DNA Binding"
cleanUpNames[grep("Cyp",cleanUpNames)] <- "CYP"
cleanUpNames[grep("Gpcr",cleanUpNames)] <- "GPCR"
endPointInfo$intended_target_family <- cleanUpNames

chemicalSummary.orig <- readRDS(file.path(pathToApp,"chemicalSummary_ACC.rds"))
chemicalSummary.orig$chnm[chemicalSummary.orig$chnm == "4-(1,1,3,3-Tetramethylbutyl)phenol"] <- "4-tert-Octylphenol"

stationINFO <- readRDS(file.path(pathToApp,"sitesOWC.rds"))

chemicalSummary <- chemicalSummary.orig %>%
  filter(endPoint %in% ep$endPoint) %>%
  data.table() %>%
  left_join(data.table(ep), by="endPoint") %>%
  data.frame() %>%
  select_("EAR","chnm","class","date","groupCol","site","endPoint","casrn") %>% 
  rename(choices = groupCol) %>%
  mutate(category=chnm) %>%
  left_join(stationINFO[,c("fullSiteID","shortName")], by=c("site"="fullSiteID")) %>%
  select(-site) %>%
  rename(site=shortName)

chemicalSummary$site[chemicalSummary$site == "Saginaw2"] <- "Saginaw"
chemicalSummary$site[chemicalSummary$site == "Cheboygan2"] <- "Cheboygan"
chemicalSummary$site[chemicalSummary$site == "Kalamazoo2"] <- "Kalamazoo"

stationINFO$shortName[stationINFO$shortName == "Kalamazoo2"] <- "Kalamazoo"
stationINFO$shortName[stationINFO$shortName == "Cheboygan2"] <- "Cheboygan"


flagsShort <- c("Borderline",  "OneAbove",
                "GainAC50", "Biochemical")
# flagsShort <- c("Borderline", "OnlyHighest", "OneAbove","Noisy",
#                 "HitCall", "GainAC50", "Biochemical")
for(i in flagsShort){
  take.out.flags <- flagDF[!flagDF[[i]],c("casn","endPoint")]
  
  chemicalSummary <- right_join(chemicalSummary, take.out.flags, 
                                by=c("casrn"="casn", "endPoint"="endPoint")) %>%
    filter(!is.na(chnm))
}

chemicalSummary <- left_join(chemicalSummary, select(flagDF, casn,endPoint, flags),
                             by=c("casrn"="casn", "endPoint"="endPoint"))

flagSum <- chemicalSummary %>%
  filter(EAR > 10^-3) %>%
  filter(!is.na(flags)) %>%
  group_by(chnm, endPoint, flags) %>%
  summarise(count = n(),
            max = max(EAR, na.rm = TRUE)) %>%
  filter(count > 10) %>%
  arrange(desc(max))

graphData <- chemicalSummary %>%
  group_by(site,date,category,class) %>%
  summarise(sumEAR=sum(EAR)) %>%
  data.frame() %>%
  group_by(site, category,class) %>%
  summarise(meanEAR=max(sumEAR)) %>%
  data.frame() 

graphData$class[graphData$class == "Human Drug, Non Prescription"] <- "Human Drug"
graphData$class[graphData$class == "Antimicrobial Disinfectant"] <- "Antimicrobial"
graphData$class[graphData$class == "Detergent Metabolites"] <- "Detergent"

orderClass <- graphData %>%
  group_by(class,category) %>%
  summarise(median = median(meanEAR[meanEAR != 0])) %>%
  data.frame() %>%
  arrange(desc(median)) %>%
  filter(!duplicated(class)) %>%
  arrange(median) 

orderChem <- graphData %>%
  group_by(category,class) %>%
  summarise(median = quantile(meanEAR[meanEAR != 0],0.5)) %>%
  data.frame() %>%
  mutate(class = factor(class, levels=orderClass$class)) %>%
  arrange(class, median)

orderClass <- mutate(orderClass, class = factor(class, levels=levels(orderChem$class)))

orderedLevels <- orderChem$category#[!is.na(orderChem$median)]


fancyNumbers2 <- function(n){
  textReturn <-  signif(n,digits = 2)
  textReturn <- as.character(textReturn)
  textReturn[length(textReturn)] <- paste(">",textReturn[length(textReturn)])
  textReturn[1] <- paste("<",textReturn[1])
  return(textReturn)
}

fancyNumbers <- function(n){
  nNoNA <- n[!is.na(n)]
  x <-gsub(pattern = "1e",replacement = "10^",x = format(nNoNA, scientific = TRUE))
  exponents <- as.numeric(sapply(strsplit(x, "\\^"), function(j) j[2]))
  base <- ifelse(exponents == 0, "1", ifelse(exponents == 1, "10","10^"))
  exponents[exponents == 0 | exponents == 1] <- ""
  textNums <- rep(NA, length(n))  
  textNums[!is.na(n)] <- paste0(base,exponents)
  
  textReturn <- parse(text=textNums)
  return(textReturn)
}

stationINFO <- stationINFO

stationINFO$Lake[stationINFO$shortName == "BlackMI"] <- "Lake Huron"
stationINFO$Lake[stationINFO$shortName == "HuronMI"] <- "Lake Erie"
stationINFO$Lake[stationINFO$shortName == "ClintonDP"] <- "Lake Erie"
stationINFO$Lake[stationINFO$shortName == "Clinton"] <- "Lake Erie"


sitesOrdered <- c("StLouis","Pigeon","Nemadji","WhiteWI","Bad","Montreal","PresqueIsle",
                  "Ontonagon","Sturgeon","Tahquamenon",
                  "Burns","IndianaHC","StJoseph","PawPaw",        
                  "Kalamazoo2","Kalamazoo","GrandMI","MilwaukeeMouth","Muskegon",      
                  "WhiteMI","Sheboygan","PereMarquette","Manitowoc",    
                  "Manistee","Fox","Oconto","Peshtigo",      
                  "Menominee","Indian","Cheboygan2","Cheboygan","Ford",         
                  "Escanaba","Manistique",
                  "ThunderBay","AuSable","Rifle",
                  "Saginaw","Saginaw2","BlackMI","Clinton","Rouge","HuronMI","Raisin","Maumee",
                  "Portage","Sandusky","HuronOH","Vermilion","BlackOH","Rocky","Cuyahoga","GrandOH",
                  "Cattaraugus","Tonawanda","Genesee","Oswego","BlackNY","Oswegatchie","Grass","Raquette","StRegis")

siteToFind <- unique(graphData$site)

siteLimits <- stationINFO %>%
  filter(shortName %in% unique(graphData$site))

siteLimits$Lake[siteLimits$Lake == "Detroit River and Lake St. Clair"] <- "Lake Michigan"
siteLimits$Lake[siteLimits$Lake == "St. Lawrence River"] <- "Lake Ontario"

siteLimits <- siteLimits %>%
  mutate(shortName = factor(shortName, levels=sitesOrdered[sitesOrdered %in% siteLimits$shortName]))


