## ----setup, include=FALSE-------------------------------------------------------------------------
options(continue=" ")
options(width=100)
library(knitr)
library(rmarkdown)
library(dplyr)
library(reshape2)


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
filePath <- file.path(packagePath, "pCodeInfo.RData")
load(file=filePath)

pCodeInfo <- pCodeInfo[pCodeInfo$casrn != "", ]

# Some chemicals have 2 rows (different units)...
# To simplify for now, we'll ignore the duplicate chemicals:
passiveData <- passiveData[!duplicated(passiveData$CAS),]

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
AC50Converted <- left_join(AC50, AC50gain)
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


## ----warning=FALSE--------------------------------------------------------------------------------
siteColumns <- grep("site",names(passiveData))
maxMinSummary <- left_join(maxMinSummary,passiveData[,-siteColumns], 
                           by=c("casn"="CAS", "Units"="Units")) 


## -------------------------------------------------------------------------------------------------
sum(maxMinSummary$MLD > maxMinSummary$minEndPoint)


## ----echo=FALSE-----------------------------------------------------------------------------------
dfToPrint <- mutate(maxMinSummary, MLD_EAR = MLD/minEndPoint) %>%
  filter(MLD_EAR > 1) %>%
  select(Chemical, minEndPoint, MLD, MLD_EAR) %>%
  arrange(desc(MLD_EAR)) %>%
  rename("Max EAR"=MLD_EAR)
  
kable(dfToPrint, digits=2)

## ----warning=FALSE--------------------------------------------------------------------------------
siteColumns <- grep("site",names(passiveData))
passiveData[,siteColumns] <- suppressWarnings(
  sapply(passiveData[,siteColumns], function(x) as.numeric(x)))

#For this analysis, we'll consider NA's to be 0 (other options exist):
passiveData[,siteColumns][is.na(passiveData[,siteColumns])] <- 0


## -------------------------------------------------------------------------------------------------
siteKey <- setNames(stationINFO$shortName, stationINFO$STAID)

dataSummary <- select(passiveData, CAS, Units) %>%
  mutate(maxMeasure = apply(passiveData[,siteColumns], 1, max, na.rm=TRUE)) %>%
  mutate(minMeasure = apply(passiveData[,siteColumns], 1, min, na.rm=TRUE)) %>%
  mutate(maxIndex = apply(passiveData[,siteColumns], 1, which.max)) %>%
  mutate(site = gsub("site","",names(passiveData[,siteColumns])[maxIndex])) %>%
  mutate(shortname = siteKey[site])

maxMinSummaryNew <- left_join(maxMinSummary, dataSummary, 
                              by=c("casn" = "CAS", "Units"="Units")) %>%
  mutate(EAR = maxMeasure/minEndPoint) %>%
  filter(EAR > 0.10) %>%
  arrange(desc(EAR)) %>%
  select(Chemical, minEndPoint, maxMeasure, EAR, shortname) %>%
  rename("Station"=shortname)
  

kable(maxMinSummaryNew, digits=2)



## ----echo=TRUE------------------------------------------------------------------------------------
endPointSummary <- select(endPoint, casn, Units) %>%  
  mutate(minEndPoint = apply(endPointData, 1, min, na.rm=TRUE)) %>%
  #Other filters would be done here (for example, only ATG, or find the closest)
  mutate(endPoint = colnames(endPointData)[apply(endPointData, 1, which.min)]) %>% 
  filter(casn %in% passiveData$CAS &
           Units %in% passiveData$Units) %>%
  slice(match(passiveData$CAS, casn)[!is.na(match(passiveData$CAS, casn))]) 
# kable(head(endPointSummary), digits=2)


## ----message=FALSE--------------------------------------------------------------------------------
commonColumns <- c("Chemical", "Units", 
                  "MLD", "MQL", "CAS", "mlWt")
oneSite <- passiveData[,c(commonColumns, "site04101500")]
oneSite <- rename(oneSite, value = site04101500)
head(oneSite)


## ----warning=FALSE--------------------------------------------------------------------------------
passiveRatio <- right_join(oneSite, endPointSummary, 
                           by=c("CAS"="casn", "Units"="Units")) %>%
  mutate(EAR = value/minEndPoint)

  

## ----echo=TRUE, eval=TRUE, warning=FALSE----------------------------------------------------------

siteColumns <- grep("site",names(passiveData))
passiveRatio <- right_join(passiveData, endPointSummary, 
                           by=c("CAS"="casn", "Units"="Units"))
passiveRatio[,siteColumns] <- passiveRatio[,siteColumns]/passiveRatio$minEndPoint

maxKey <- setNames(apply(passiveRatio[,siteColumns], 2, max),
                   gsub("site","", names(passiveRatio[,siteColumns])))
chemKey <- setNames(passiveRatio$Chemical[apply(passiveRatio[,siteColumns], 2, which.max)],
                   gsub("site","", names(passiveRatio[,siteColumns])))

endKey <- setNames(passiveRatio$endPoint[apply(passiveRatio[,siteColumns], 2, which.max)],
                   gsub("site","", names(passiveRatio[,siteColumns])))


maxRatioBySite <-  data.frame(site=names(passiveData)[siteColumns],
                             stringsAsFactors=FALSE) %>%
  mutate(shortName = siteKey[gsub("site","",site)]) %>%
  mutate(maxRatio = maxKey[gsub("site","",site)]) %>%
  mutate(Chemical = chemKey[gsub("site","",site)]) %>%
  mutate(Endpoint = endKey[gsub("site","",site)]) %>%
  filter(maxRatio > 50) %>%
  arrange(desc(maxRatio)) %>%
  select(shortName, maxRatio, Chemical, Endpoint) %>%
  rename(Station=shortName, "Max Ratio"=maxRatio, "End Point"=Endpoint)
 
kable(maxRatioBySite, digits=2, row.names = FALSE)


## ----warning=FALSE, echo=FALSE, results='asis'----------------------------------------------------
infoColumns <- c("Chemical", "CAS")
endpointNames <- names(endPoint)
endpointNames <- endpointNames[!(endpointNames %in% c("casn","Units","mlWt","conversion"))]

for(i in siteColumns){
  oneSite <- passiveData[,infoColumns]
  oneSite$value <- passiveData[,i]
  
  oneSiteLong <- filter(oneSite, value != 0) %>%
    rename(measuredValue = value) %>%
    left_join(pCodeInfo[c("parameter_cd","casrn","class","parameter_nm")], by=c("CAS"="casrn"))%>%
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
    arrange(chnm, desc(EAR)) %>%
    rename(type=intended_target_type, type_sub=intended_target_type_sub,
           family=intended_target_family, family_sub=intended_target_family_sub) %>%
    select(-contains("intended_target_"))
  
  oneSiteLong <- unique(oneSiteLong)


  if(nrow(oneSiteLong) != 0){
    cat("\n###",siteKey[gsub("site","",names(passiveData)[i])])
    
    print(kable(oneSiteLong, digits=2, caption = siteKey[gsub("site","",names(passiveData)[i])], row.names = FALSE))
  }
  
  
}


## ----message=FALSE, results='asis'----------------------------------------------------------------
infoColumns <- c("Chemical", "CAS", "Units", "MLD", "MQL", "mlWt")


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
  group_by(Chemical, endPoint) %>%
  summarize(minEAR=min(EAR), maxEAR=max(EAR), hits=length(EAR), nSites=ncol(passiveData)) %>%
  left_join(endPointInfo, by=c("endPoint"="assay_component_endpoint_name")) %>%
  select(Chemical, minEAR, maxEAR, hits, nSites, endPoint, contains("intended_target_")) %>%
  arrange(Chemical, desc(maxEAR)) %>%
  rename(type=intended_target_type, type_sub=intended_target_type_sub,
         family=intended_target_family, family_sub=intended_target_family_sub) %>%
  select(-contains("intended_target_")) %>%
  arrange(desc(maxEAR))
  
  print(kable(chemicalSummary, digits=2, caption = "Passive chemical summary", row.names = FALSE))


## ----echo=FALSE-----------------------------------------------------------------------------------

#Clear out local enviornment:
rm(AC50, dataSummary, dfToPrint,endPoint,
   endPointData, maxMinSummary, maxMinSummaryNew, maxRatioBySite,
   oneSite, passiveData,  commonColumns, filePath,
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


## ----warning=FALSE--------------------------------------------------------------------------------
detLevels  <- data.frame(dl_pcode=names(detectionLimits),
                         minDetLevel=apply(detectionLimits, 2, min),
                         row.names=NULL, stringsAsFactors=FALSE) %>%
  mutate(pCode=sapply(strsplit(dl_pcode, "_"), function(x) x[2])) %>%
  left_join(pCodeInfo[c("parameter_cd","casrn")], by=c("pCode"="parameter_cd")) %>%
  select(casrn, pCode, minDetLevel) %>%
  right_join(maxMinSummary, by=c("casrn"="casn")) %>%
  mutate(EAR=minDetLevel/minEndPoint) %>%
  arrange(desc(EAR)) %>%
  filter(EAR > 0.1)

kable(detLevels, digits = 3)


## -------------------------------------------------------------------------------------------------
qualColumns <- grep("qualifier", names(waterSamples))
valColumns <- grep("valueToUse", names(waterSamples))

waterData <- waterSamples[,valColumns]
waterData[waterSamples[,qualColumns] == "<"] <- NA


## ----echo=FALSE, warning=FALSE--------------------------------------------------------------------

dataSummary  <- suppressWarnings(data.frame(pcode=names(waterData),
                         maxValue=apply(waterData, 2, max, na.rm=TRUE, na.action=NA),
                         row.names=NULL, stringsAsFactors=FALSE) ) %>%
  mutate(pCode=sapply(strsplit(pcode, "_"), function(x) x[2])) %>%
  left_join(pCodeInfo[c("parameter_cd","casrn")], by=c("pCode"="parameter_cd")) %>%
  select(casrn, pCode, maxValue) %>%
  right_join(maxMinSummary, by=c("casrn"="casn")) %>%
  mutate(EAR=maxValue/minEndPoint) %>%
  arrange(desc(EAR)) %>%
  filter(EAR > 0.1)
  

kable(dataSummary, digits=3)


## ----message=FALSE, results='asis'----------------------------------------------------------------

for(i in unique(waterSamples$site)){
  
  oneSite <- cbind(waterSamples[,1:2],waterData) %>%
    filter(site==i) 
  
  valColumns <- grep("valueToUse",names(oneSite))
  shortName <- stationINFO$shortName[i == stationINFO$fullSiteID]
  
  
  
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
    arrange(desc(maxEAR))

  if(nrow(oneSiteLong2) > 0){
    cat("\n\n###", shortName, "\n")
    print(kable(oneSiteLong2, digits=3,caption = shortName, row.names = FALSE))
  }
  

}


## ----message=FALSE, results='asis'----------------------------------------------------------------
valColumns <- grep("valueToUse",names(oneSite))

chemicalSummary <- cbind(waterSamples[,1:2],waterData) %>%
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
  group_by(chnm, endPoint) %>%
  summarize(minEAR=min(EAR), maxEAR=max(EAR), hits=length(EAR), nSites=length(unique(site))) %>%
  left_join(endPointInfo, by=c("endPoint"="assay_component_endpoint_name")) %>%
  select(chnm, minEAR, maxEAR, hits, nSites, endPoint, contains("intended_target_")) %>%
  arrange(chnm, desc(maxEAR)) %>%
  rename(type=intended_target_type, type_sub=intended_target_type_sub,
         family=intended_target_family, family_sub=intended_target_family_sub) %>%
  select(-contains("intended_target_")) %>%
  arrange(desc(maxEAR))
  
  print(kable(chemicalSummary, digits=2, caption = "Water Sample Chemical Summary", row.names = FALSE))


## ----fig.width=7, fig.height=10, warning=FALSE, echo=FALSE, fig.cap="Colors represent different endpoints, dots represent measured water sample values"----
library(ggplot2)

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
  filter(endPointValue <= 50)

# ggplot(ep[1:50,], aes(variable)) + geom_bar() + facet_wrap(~ chnm)
# ggplot(ep[1:50,], aes(variable, fill=chnm)) + geom_bar() # + scale_fill_brewer()

ggplot(data, aes(x = chnm, y = endPointValue, fill=endPoint)) +
    geom_bar(stat='identity',guide=FALSE,colour="grey") + 
    ylim(0,10) +   
    theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") +
    labs(x = "Chemical", y="Concentration") + 
    geom_point(aes(x=chnm, y=measuredValue))


