## ----setup, include=FALSE-------------------------------------------------------------------------
options(continue=" ")
options(width=100)
library(knitr)
library(rmarkdown)
library(dplyr)
library(reshape2)
library(DT)


## -------------------------------------------------------------------------------------------------

library("webchem")

molweight <- cir_query('3380-34-5', "mw")
molweight


## ----echo=FALSE, eval=TRUE------------------------------------------------------------------------

packagePath <- system.file("extdata", package="toxEval")
filePath <- file.path(packagePath, "passiveData.RData")
load(file=filePath)
#Site Info:
filePath <- file.path(packagePath, "stationINFO.RData")
load(file=filePath)
#parameter info:
filePath <- file.path(packagePath, "pCodeInfo.RData")
load(file=filePath)

#passive class info:
passiveCAS <- readRDS(file.path(packagePath, "passiveCAS.rds"))

pCodeInfo <- pCodeInfo[pCodeInfo$casrn != "", ]

# Some chemicals have 2 rows (different units)...
# To simplify for now, we'll ignore the duplicate chemicals:


passiveData <- left_join(passiveData, passiveCAS, by=c("CAS", "Chemical"))
passiveData <- passiveData[!duplicated(passiveData$CAS),]
passiveData$class[passiveData$Chemical == "Aspirin"] <-"Pharmaceuticals"

passiveData <- passiveData[,-which(names(passiveData) == "site04249000.1")]

head(passiveData[,1:7])

# Unique units:
unique(passiveData$Units)


## ----message=FALSE--------------------------------------------------------------------------------
library(toxEval)

# AC50 data provided in the toxEval package:
AC50gain <- AC50gain
head(AC50gain[,1:5])


## ----message=FALSE--------------------------------------------------------------------------------
library(dplyr)

unitConversion <- setNames(c(10^6, 10^3), c("pg/L", "ng/L") )

AC50 <- right_join(AC50gain[,c("casn"), drop=FALSE], 
                   passiveData[,c("CAS", "Units", "mlWt")],
                   by= c("casn" = "CAS")) %>%
  rename(desiredUnits = Units) %>%
  filter(!is.na(mlWt)) %>%
  mutate(conversion = unitConversion[desiredUnits] * mlWt) %>%
  select(casn, desiredUnits, mlWt, conversion)


## ----echo=FALSE-----------------------------------------------------------------------------------
kable(head(AC50), digits=2, row.names = FALSE)


## -------------------------------------------------------------------------------------------------
AC50Converted <- left_join(AC50, AC50gain, by="casn")
infoColumns <- c("casn", "chnm", "desiredUnits","mlWt", "conversion", "code","chid")

endPointData <- AC50Converted[,!(names(AC50Converted) %in% infoColumns)]
endPointData <- 10^endPointData
endPointData <- endPointData * AC50Converted$conversion
endPoint <- cbind(AC50, data.frame(endPointData)) %>%
  rename(Units=desiredUnits) %>%
  filter(rowSums(is.na(endPointData)) != ncol(endPointData))

endPointData <- endPointData[rowSums(is.na(endPointData)) != 
                               ncol(endPointData),]

## ----warning=FALSE--------------------------------------------------------------------------------
maxEndPoints <- apply(endPointData, 1, max, na.rm=TRUE) 
minEndPoints <- apply(endPointData, 1, min, na.rm=TRUE)

maxMinSummary <- cbind(endPoint[,c("casn", "Units")], 
                       maxEndPoint=maxEndPoints, 
                       minEndPoint=minEndPoints) %>%
  arrange(desc(casn))


## ----warning=FALSE, echo=FALSE--------------------------------------------------------------------
siteColumns <- grep("site",names(passiveData))
maxMinSummary <- left_join(maxMinSummary,passiveData[,-siteColumns], 
                           by=c("casn"="CAS", "Units"="Units")) 


## ----echo=FALSE-----------------------------------------------------------------------------------
sum(maxMinSummary$MLD > maxMinSummary$minEndPoint)


## ----echo=FALSE-----------------------------------------------------------------------------------
dfToPrint <- mutate(maxMinSummary, MLD_EAR = MLD/minEndPoint) %>%
  filter(MLD_EAR > 0.01) %>%
  select(Chemical, minEndPoint, MLD, MLD_EAR, class) %>%
  arrange(desc(MLD_EAR)) %>%
  rename("Max EAR"=MLD_EAR)
  
datatable(dfToPrint, rownames=FALSE) %>% 
    formatRound(c("Max EAR","MLD","minEndPoint"), digits = 3)


## ----warning=FALSE--------------------------------------------------------------------------------
siteColumns <- grep("site",names(passiveData))
passiveData[,siteColumns] <- suppressWarnings(
  sapply(passiveData[,siteColumns], function(x) as.numeric(x)))

#For this analysis, we'll consider NA's to be 0 (other options exist):
passiveData[,siteColumns][is.na(passiveData[,siteColumns])] <- 0


## ----message=FALSE, results='asis', echo=FALSE----------------------------------------------------
infoColumns <- c("Chemical", "CAS", "Units", "MLD", "MQL", "mlWt","class")

chemicalSummary <- melt(passiveData, id.vars = infoColumns) %>%
  mutate(variable=as.character(variable)) %>%
  filter(value != 0) %>%
  rename(measuredValue=value, site=variable) %>%
  right_join(endPoint, by=c("CAS"="casn", "Units"="Units", "mlWt"="mlWt")) %>%
  melt(id.vars = c(infoColumns, "site", "measuredValue","conversion")) %>% 
  mutate(variable=as.character(variable)) %>%
  rename(endPointValue=value, endPoint=variable) %>%
  filter(!is.na(endPointValue)) %>%
  mutate(EAR=measuredValue/endPointValue) %>%
  filter(EAR > 0.1) %>%
  group_by(Chemical, endPoint, class) %>%
  summarize(minEAR=min(EAR), maxEAR=max(EAR), hits=length(EAR), nSites=ncol(passiveData)) %>%
  left_join(endPointInfo, by=c("endPoint"="assay_component_endpoint_name")) %>%
  select(Chemical, minEAR, maxEAR, hits, nSites, endPoint, class, contains("intended_target_")) %>%
  arrange(Chemical, desc(maxEAR)) %>%
  rename(type=intended_target_type, type_sub=intended_target_type_sub,
         family=intended_target_family, family_sub=intended_target_family_sub) %>%
  select(-contains("intended_target_")) %>% 
    data.frame %>%
  arrange(desc(maxEAR)) 

  chemicalSummary <- unique(chemicalSummary) #some endpoints have multiple genes....so there are duplicated rows, there wouldn't be if we didn't have select(-contains("intended_target_"))

  datatable(chemicalSummary, rownames = FALSE, 
            options = list(pageLength = 10)) %>% 
    formatRound(c("minEAR","maxEAR"), digits = 2)


## ----message=FALSE, results='asis', echo=FALSE----------------------------------------------------
infoColumns <- c("Chemical", "CAS", "Units", "MLD", "MQL", "mlWt","class")
siteKey <- setNames(stationINFO$shortName, stationINFO$STAID)

siteSummary <- melt(passiveData, id.vars = infoColumns) %>%
  mutate(variable=as.character(variable)) %>%
  filter(value != 0) %>%
  rename(measuredValue=value, site=variable) %>%
  right_join(endPoint, by=c("CAS"="casn", "Units"="Units", "mlWt"="mlWt")) %>%
  melt(id.vars = c(infoColumns, "site", "measuredValue","conversion")) %>% 
  mutate(variable=as.character(variable)) %>%
  rename(endPointValue=value, endPoint=variable) %>%
  filter(!is.na(endPointValue)) %>%
  mutate(EAR=measuredValue/endPointValue) %>%
  filter(EAR > 0.1) %>%
  mutate(site=siteKey[gsub("site","",site)])%>%
  group_by(site, endPoint) %>%
  summarize(minEAR=min(EAR), maxEAR=max(EAR), hits=length(EAR), nChem=length(Chemical), nEndPoints=length(unique(endPoint))) %>%
  select(site, minEAR, maxEAR, hits, nChem, nEndPoints)  %>%
    data.frame %>%
  arrange(desc(maxEAR))  

  datatable(siteSummary, rownames = FALSE, 
            options = list(pageLength = 10)) %>% 
    formatRound(c("minEAR","maxEAR"), digits = 2)


## ----warning=FALSE, echo=FALSE, results='asis'----------------------------------------------------
infoColumns <- c("Chemical", "CAS", "Units", "MLD", "MQL", "mlWt","class")

endpointNames <- names(endPoint)
endpointNames <- endpointNames[!(endpointNames %in% c("casn","Units","mlWt","conversion"))]

reversesiteKey <- setNames(paste0("site",stationINFO$STAID), stationINFO$shortName)

for(i in siteSummary$site){
  oneSite <- passiveData[,infoColumns]
  index <- reversesiteKey[i]
  oneSite$value <- passiveData[,index]
  
  oneSiteLong <- filter(oneSite, value != 0) %>%
    rename(measuredValue = value) %>%
#     left_join(pCodeInfo[c("parameter_cd","casrn","parameter_nm")], by=c("CAS"="casrn")) %>%
    select(Chemical, CAS, class, measuredValue) %>%
    right_join(endPoint, by=c("CAS"="casn")) %>%
    select(-mlWt, -conversion, -CAS,  -Units) %>%
    rename(chnm=Chemical) %>%
    melt(id.vars = c("measuredValue", "chnm", "class")) %>%
    mutate(variable=as.character(variable)) %>%
    rename(endPointValue=value, endPoint=variable) %>%
    filter(!is.na(endPointValue)) %>%
    mutate(EAR=measuredValue/endPointValue) %>%
    filter(EAR > 0.1) %>%
    left_join(endPointInfo, by=c("endPoint"="assay_component_endpoint_name")) %>%
    select(chnm, EAR, class, endPoint, contains("intended_target_")) %>%
    rename(type=intended_target_type, type_sub=intended_target_type_sub,
           family=intended_target_family, family_sub=intended_target_family_sub) %>%
    select(-contains("intended_target_")) %>%
    arrange(desc(EAR))
  
  oneSiteLong <- unique(oneSiteLong)
    

  if(nrow(oneSiteLong) != 0){
    cat("\n###",i)
    
#     datatable(oneSiteLong, rownames = FALSE) %>% formatRound(c("EAR"), digits = 2)
    print(kable(oneSiteLong, digits=2, caption = i, row.names = FALSE))
  }
  
  
}


## ----fig.width=10, fig.height=10, warning=FALSE, echo=FALSE, fig.cap="Colors represent different endpoints, dots represent measured passive water sample values"----
library(ggplot2)
infoColumns <- c("Chemical", "CAS", "Units", "MLD", "MQL", "mlWt","class")

maxY <- 50000
dataPassive <- melt(passiveData, id.vars = infoColumns) %>%
  mutate(variable=as.character(variable)) %>%
  filter(value != 0) %>%
  rename(measuredValue=value, site=variable) %>%
  right_join(endPoint, by=c("CAS"="casn", "Units"="Units", "mlWt"="mlWt")) %>%
  melt(id.vars = c(infoColumns, "site", "measuredValue","conversion")) %>% 
  mutate(variable=as.character(variable)) %>%
  rename(endPointValue=value, endPoint=variable) %>%
  filter(!is.na(endPointValue)) %>%
  mutate(endPointValue=ifelse(endPointValue>maxY, maxY, endPointValue))

# ggplot(ep[1:50,], aes(variable)) + geom_bar() + facet_wrap(~ chnm)
# ggplot(ep[1:50,], aes(variable, fill=chnm)) + geom_bar() # + scale_fill_brewer()

ggplot(dataPassive, aes(x = Chemical, y = endPointValue, fill=endPoint)) +
    geom_bar(stat='identity',guide=FALSE,colour="grey") + 
    ylim(0,maxY) +   
    theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.25), legend.position = "none") +
    labs(x = "Chemical", y="Concentration [ng/L]") + 
    geom_point(aes(x=Chemical, y=measuredValue))


## ----echo=FALSE-----------------------------------------------------------------------------------

#Clear out local enviornment:
rm(AC50, dataSummary, dfToPrint,endPoint,
   endPointData, maxMinSummary, maxMinSummaryNew,
   oneSite, passiveData,  filePath,
   infoColumns, maxEndPoints,
   minEndPoints, packagePath, siteColumns,unitConversion)


## ----echo=TRUE, eval=TRUE-------------------------------------------------------------------------
packagePath <- system.file("extdata", package="toxEval")
filePath <- file.path(packagePath, "waterSamples.RData")
load(file=filePath)

head(waterSamples[,1:10])


## ----message=FALSE, eval=FALSE--------------------------------------------------------------------
#  library(dataRetrieval)
#  
#  waterSamplePCodes <- names(waterSamples)[grep("valueToUse", names(waterSamples))]
#  waterSamplePCodes <- sapply(strsplit(waterSamplePCodes, "_"), function(x) x[2])
#  
#  pCodeInfo <- readNWISpCode(waterSamplePCodes)
#  
#  unique(pCodeInfo$parameter_units)
#  
#  library(webchem)
#  pCodeInfo$mlWt <- rep(NA, nrow(pCodeInfo))
#  pCodeInfo$mlWt[pCodeInfo$casrn != ""] <- sapply(pCodeInfo$casrn[pCodeInfo$casrn != ""],
#                  function(x) cir_query(x, "mw", first = TRUE))
#  pCodeInfo$mlWt <- as.numeric(pCodeInfo$mlWt)
#  

## ----message=FALSE, echo=FALSE--------------------------------------------------------------------

waterSamplePCodes <- names(waterSamples)[grep("valueToUse", names(waterSamples))]
waterSamplePCodes <- sapply(strsplit(waterSamplePCodes, "_"), function(x) x[2])




## -------------------------------------------------------------------------------------------------

AC50 <- left_join(AC50gain, pCodeInfo[,c("casrn", "parameter_units", "mlWt")],
                   by= c("casn"="casrn")) %>%
  filter(!is.na(parameter_units)) %>%
  rename(desiredUnits = parameter_units) %>%
  mutate(conversion = mlWt) %>%
  select(casn, chnm, desiredUnits, mlWt, conversion)


## -------------------------------------------------------------------------------------------------
AC50Converted <- left_join(AC50, AC50gain)
infoColumns <- c("casn", "chnm", "desiredUnits","mlWt", "conversion", "code","chid")

endPointData <- AC50Converted[,!(names(AC50Converted) %in% infoColumns)]
endPointData <- 10^endPointData
endPointData <- endPointData * AC50Converted$conversion
endPoint <- cbind(AC50, data.frame(endPointData))
endPoint <- rename(endPoint, Units=desiredUnits)


## ----warning=FALSE--------------------------------------------------------------------------------
maxEndPoints <- apply(endPointData, 1, max, na.rm=TRUE) 
minEndPoints <- apply(endPointData, 1, min, na.rm=TRUE)

maxMinSummary <- cbind(endPoint[,c("chnm", "casn", "Units")], 
                       maxEndPoint=maxEndPoints, 
                       minEndPoint=minEndPoints)


## -------------------------------------------------------------------------------------------------
pCodeSummary <- select(pCodeInfo, srsname, casrn, parameter_units, mlWt)

detectionLimits <- waterSamples[,grep("detectionLimit", names(waterSamples))]
detectionLimits <- na.omit(detectionLimits)

# Check if any change over time:
all(apply(detectionLimits, 2, function(x) length(unique(x))==1))

#They do!


## ----warning=FALSE, echo=FALSE--------------------------------------------------------------------
detLevels  <- data.frame(dl_pcode=names(detectionLimits),
                         minDetLevel=apply(detectionLimits, 2, min),
                         row.names=NULL, stringsAsFactors=FALSE) %>%
  mutate(pCode=sapply(strsplit(dl_pcode, "_"), function(x) x[2])) %>%
  left_join(pCodeInfo[c("parameter_cd","casrn","class")], by=c("pCode"="parameter_cd")) %>%
  select(casrn, class, minDetLevel) %>%
  right_join(maxMinSummary, by=c("casrn"="casn")) %>%
  mutate(EAR=minDetLevel/minEndPoint) %>%
  arrange(desc(EAR)) %>%
  filter(EAR > 0.1) %>%
  rename(Chemical=chnm) %>%
  select(Chemical, EAR, minDetLevel, minEndPoint, class)

datatable(detLevels, rownames = FALSE)  %>% 
    formatRound(c("EAR","minDetLevel", "minEndPoint"), digits = 3)


## -------------------------------------------------------------------------------------------------
qualColumns <- grep("qualifier", names(waterSamples))
valColumns <- grep("valueToUse", names(waterSamples))

waterData <- waterSamples[,valColumns]
waterData[waterSamples[,qualColumns] == "<"] <- NA


## ----message=FALSE, results='asis', echo=FALSE----------------------------------------------------
newSiteKey <- setNames(stationINFO$shortName, stationINFO$fullSiteID)

wData <- cbind(waterSamples[,1:2],waterData)
valColumns <- grep("valueToUse",names(wData))

chemicalSummary <- wData %>%
  melt(id.vars=c("ActivityStartDateGiven","site")) %>%
  mutate(variable=as.character(variable)) %>%
  filter(!is.na(value)) %>%
  rename(measuredValue=value, pCode=variable) %>%
  mutate(pCode=sapply(strsplit(pCode, "_"), function(x) x[2])) %>%
  left_join(pCodeInfo[c("parameter_cd","casrn","class","srsname")], by=c("pCode"="parameter_cd")) %>%
  select(srsname,casrn, class, measuredValue, site) %>%
  right_join(endPoint, by=c("casrn"="casn")) %>%
  select(-mlWt, -conversion, -casrn,  -Units, -srsname) %>%
  melt(id.vars = c("class", "site", "measuredValue","chnm")) %>% 
  filter(!is.na(value)) %>%
  mutate(variable=as.character(variable)) %>%
  rename(endPointValue=value, endPoint=variable) %>%
  mutate(EAR=measuredValue/endPointValue) %>%
  filter(EAR > 0.1) %>%
  group_by(chnm, endPoint, class) %>%
  summarize(minEAR=min(EAR), maxEAR=max(EAR), hits=length(EAR), nSites=length(unique(site))) %>%
  left_join(endPointInfo, by=c("endPoint"="assay_component_endpoint_name")) %>%
  select(chnm, minEAR, maxEAR, hits, nSites, endPoint, class, contains("intended_target_")) %>%
  arrange(chnm, desc(maxEAR)) %>%
  rename(type=intended_target_type, type_sub=intended_target_type_sub,
         family=intended_target_family, family_sub=intended_target_family_sub) %>%
  select(-contains("intended_target_"))  %>%
  data.frame %>%
  arrange(desc(maxEAR))
  

  datatable(chemicalSummary, rownames = FALSE, 
            options = list(pageLength = 10)) %>% 
    formatRound(c("minEAR","maxEAR"), digits = 2)
#   print(kable(chemicalSummary, digits=2, caption = "Water Sample Chemical Summary", row.names = FALSE))


## ----message=FALSE, results='asis', echo=FALSE----------------------------------------------------

siteSummary <- wData %>%
  melt(id.vars=c("ActivityStartDateGiven","site")) %>%
  mutate(variable=as.character(variable)) %>%
  filter(!is.na(value)) %>%
  rename(measuredValue=value, pCode=variable) %>%
  mutate(pCode=sapply(strsplit(pCode, "_"), function(x) x[2])) %>%
  left_join(pCodeInfo[c("parameter_cd","casrn","class","srsname")], by=c("pCode"="parameter_cd")) %>%
  select(srsname,casrn, class, measuredValue, site) %>%
  right_join(endPoint, by=c("casrn"="casn")) %>%
  select(-mlWt, -conversion, -casrn,  -Units, -srsname) %>%
  melt(id.vars = c("class", "site", "measuredValue","chnm")) %>% 
  filter(!is.na(value)) %>%
  mutate(variable=as.character(variable)) %>%
  rename(endPointValue=value, endPoint=variable) %>%
  mutate(EAR=measuredValue/endPointValue) %>%
  filter(EAR > 0.1) %>%
  mutate(site=newSiteKey[site]) %>%
  group_by(site) %>%
  summarize(minEAR=min(EAR), maxEAR=max(EAR), hits=length(EAR), nChem=length(unique(chnm)), nEndPoints=length(unique(endPoint))) %>%
  select(site, minEAR, maxEAR, hits, nChem, nEndPoints) %>%
  data.frame %>%
  arrange(desc(hits))
  

  datatable(siteSummary, rownames = FALSE, 
            options = list(pageLength = 10)) %>% 
    formatRound(c("minEAR","maxEAR"), digits = 2)


## ----message=FALSE, results='asis', echo=FALSE----------------------------------------------------

for(i in siteSummary$site){
  
  oneSite <- wData %>%
    mutate(site = newSiteKey[site]) %>%
    filter(site == i) 
  
  valColumns <- grep("valueToUse",names(oneSite))
  
  oneSiteLong2 <- melt(oneSite[,valColumns]) %>%
    mutate(variable=as.character(variable)) %>%
    mutate(pCode=sapply(strsplit(variable, "_"), function(x) x[2])) %>%
    select(-variable) %>%
    filter(!is.na(value)) %>%
    rename(measuredValue=value) %>%
    left_join(pCodeInfo[c("parameter_cd","casrn","class")], by=c("pCode"="parameter_cd")) %>%
    select(casrn, class, measuredValue) %>%
    right_join(endPoint, by=c("casrn"="casn")) %>%
    select(-mlWt, -conversion, -casrn,  -Units) %>%
    melt(id.vars = c("measuredValue", "chnm", "class")) %>%
    mutate(variable=as.character(variable)) %>%
    rename(endPointValue=value, endPoint=variable) %>%
    filter(!is.na(endPointValue)) %>%
    mutate(EAR=measuredValue/endPointValue) %>%
    filter(EAR > 0.1) %>%
    left_join(endPointInfo, by=c("endPoint"="assay_component_endpoint_name")) %>%
    select(chnm, EAR, class, endPoint, contains("intended_target_")) %>%
    arrange(chnm, desc(EAR)) %>%
    group_by(chnm, endPoint, class, intended_target_type, intended_target_type_sub, 
             intended_target_family, intended_target_family_sub) %>%
    rename(type=intended_target_type, type_sub=intended_target_type_sub,
           family=intended_target_family, family_sub=intended_target_family_sub) %>%
    summarize(minEAR=min(EAR), maxEAR=max(EAR), hits=length(EAR), count=nrow(oneSite)) %>%
    select(chnm, minEAR, maxEAR, hits, count, endPoint, class, type, type_sub, family, family_sub) %>%
    data.frame %>%
    arrange(desc(maxEAR))

  if(nrow(oneSiteLong2) > 0){
    cat("\n\n###", i, "\n")
    print(kable(oneSiteLong2, digits=3,caption = i, row.names = FALSE))
  }
  

}


## ----fig.width=10, fig.height=10, warning=FALSE, echo=FALSE, fig.cap="Colors represent different endpoints, dots represent measured water sample values"----

maxY <- 10
data <- cbind(waterSamples[,1:2],waterData) %>%
  melt(id.vars=c("ActivityStartDateGiven","site")) %>%
  mutate(variable=as.character(variable)) %>%
  filter(!is.na(value)) %>%
  rename(measuredValue=value, pCode=variable) %>%
  mutate(pCode=sapply(strsplit(pCode, "_"), function(x) x[2])) %>%
  left_join(pCodeInfo[c("parameter_cd","casrn","class","srsname")], by=c("pCode"="parameter_cd")) %>%
  select(srsname,casrn, class, measuredValue, site) %>%
  right_join(endPoint, by=c("casrn"="casn")) %>%
  select(-mlWt, -conversion, -casrn,  -Units, -srsname) %>%
  melt(id.vars = c("class", "site", "measuredValue","chnm")) %>% 
  filter(!is.na(value)) %>%
  mutate(variable=as.character(variable)) %>%
  rename(endPointValue=value, endPoint=variable) %>%
  mutate(endPointValue=ifelse(endPointValue <= maxY, maxY, endPointValue)) 

# ggplot(ep[1:50,], aes(variable)) + geom_bar() + facet_wrap(~ chnm)
# ggplot(ep[1:50,], aes(variable, fill=chnm)) + geom_bar() # + scale_fill_brewer()

ggplot(data, aes(x = chnm, y = endPointValue, fill=endPoint)) +
    geom_bar(stat='identity',guide=FALSE,colour="grey") + 
    ylim(0,10) +   
    theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.25), legend.position = "none") +
    labs(x = "Chemical", y="Concentration") + 
    geom_point(aes(x=chnm, y=measuredValue))




## ----fig.width=10, fig.height=7, warning=FALSE, echo=FALSE, fig.cap="Boxplot of all EARS over 0.1"----

data2 <- cbind(waterSamples[,1:2],waterData) %>%
  melt(id.vars=c("ActivityStartDateGiven","site")) %>%
  mutate(variable=as.character(variable)) %>%
  filter(!is.na(value)) %>%
  rename(measuredValue=value, pCode=variable) %>%
  mutate(pCode=sapply(strsplit(pCode, "_"), function(x) x[2])) %>%
  left_join(pCodeInfo[c("parameter_cd","casrn")], by=c("pCode"="parameter_cd")) %>%
  select(casrn, measuredValue, site) %>%
  right_join(endPoint, by=c("casrn"="casn")) %>%
  select(-mlWt, -conversion, -casrn,  -Units) %>%
  melt(id.vars = c("site", "measuredValue","chnm")) %>% 
  filter(!is.na(value)) %>%
  mutate(variable=as.character(variable)) %>%
  rename(endPointValue=value, endPoint=variable) %>%
  mutate(EAR = measuredValue/endPointValue) %>%
  filter(EAR >= 0.1)

ggplot(data2, aes(x = chnm, y = EAR)) +
    geom_boxplot() + 
    ylim(0,3) +    
    theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.25), legend.position = "none") +
    labs(x = "Chemical", y="EAR")



## ----fig.width=10, fig.height=7, warning=FALSE, echo=FALSE, fig.cap="Boxplot of all max EARS per site per chemical (filtered over 0.1)"----

data3 <- group_by(data2, site, chnm) %>%
  mutate(maxEAR = max(EAR))

ggplot(data3, aes(x = chnm, y = maxEAR)) +
    geom_boxplot() + 
    ylim(0,3) +    
    theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.25)) +
    labs(x = "Chemical", y="maxEAR")+ 
    geom_point(aes(x=chnm, y=EAR, color="red"))


## ----fig.width=10, fig.height=7, warning=FALSE, echo=FALSE, fig.cap="Boxplot of  EARS per site (filtered over 0.1, max per endpoint)"----

data4 <- mutate(data2, site = newSiteKey[site]) %>%
  group_by(site, chnm, endPoint) %>%
  mutate(maxEAR = max(EAR))
  

ggplot(data4, aes(x = site, y = maxEAR)) +
    geom_boxplot() + 
    ylim(0,3) +    
    theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.25), legend.position = "none") +
    labs(x = "Site", y="maxEAR")+ 
    geom_point(aes(x=site, y=EAR, color="red") )


