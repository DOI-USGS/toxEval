## ----setup, include=FALSE-------------------------------------------------------------------------
options(continue=" ")
options(width=100)
library(knitr)
library(rmarkdown)



## -------------------------------------------------------------------------------------------------

library("webchem")

molweight <- cir_query('3380-34-5', "mw")
molweight


## ----echo=TRUE, eval=TRUE-------------------------------------------------------------------------

packagePath <- system.file("extdata", package="toxEval")
filePath <- file.path(packagePath, "passiveData.RData")
load(file=filePath)

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

AC50 <- left_join(AC50gain, passiveData[,c("CAS", "Units", "mlWt")],
                   by= c("casn" = "CAS")) %>%
  filter(!is.na(Units)) %>%
  rename(desiredUnits = Units) %>%
  mutate(conversion = unitConversion[desiredUnits] * mlWt) %>%
  select(casn, chnm, desiredUnits, mlWt, conversion)


## ----echo=FALSE-----------------------------------------------------------------------------------
kable(head(AC50), digits=2, row.names = FALSE)


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


## ----warning=FALSE--------------------------------------------------------------------------------
siteColumns <- grep("site",names(passiveData))
maxMinSummary <- left_join(maxMinSummary,passiveData[,-siteColumns], 
                           by=c("casn"="CAS", "Units"="Units"))


## -------------------------------------------------------------------------------------------------
sum(maxMinSummary$MLD > maxMinSummary$minEndPoint)


## ----echo=FALSE-----------------------------------------------------------------------------------
dfToPrint <- mutate(maxMinSummary, MLD_EAR = MLD/minEndPoint) %>%
  filter(MLD_EAR > 0.01) %>%
  select(chnm, minEndPoint, MLD, MLD_EAR) %>%
  arrange(desc(MLD_EAR)) %>%
  rename("Maximum Ratio"=MLD_EAR)
  
kable(dfToPrint, digits=2)

## ----warning=FALSE--------------------------------------------------------------------------------
siteColumns <- grep("site",names(passiveData))
passiveData[,siteColumns] <- sapply(passiveData[,siteColumns], function(x) as.numeric(x))

#For this analysis, we'll consider NA's to be 0 (other options exist):
passiveData[,siteColumns][is.na(passiveData[,siteColumns])] <- 0


## ----echo=FALSE-----------------------------------------------------------------------------------
dataSummary <- select(passiveData, CAS, Units) %>%
  mutate(maxMeasure = apply(passiveData[,siteColumns], 1, max, na.rm=TRUE)) %>%
  mutate(minMeasure = apply(passiveData[,siteColumns], 1, min, na.rm=TRUE)) 

maxMinSummaryNew <- left_join(maxMinSummary, dataSummary, 
                              by=c("casn" = "CAS", "Units"="Units")) %>%
  mutate(EAR = maxMeasure/minEndPoint) %>%
  filter(EAR > 0.01) %>%
  arrange(desc(EAR)) %>%
  select(Chemical, minEndPoint, maxMeasure, EAR) 
  

kable(maxMinSummaryNew, digits=2)



## ----message=FALSE--------------------------------------------------------------------------------
commonColumns <- c("Chemical", "Units", 
                  "MLD", "MQL", "CAS", "mlWt")
oneSite <- passiveData[,c(commonColumns, "site04101500")]
oneSite <- rename(oneSite, value = site04101500)
head(oneSite)


## -------------------------------------------------------------------------------------------------
casnRow <- setNames(1:nrow(oneSite),oneSite$CAS)
ratioPassive <- oneSite[casnRow[endPoint$casn],"value"] / endPointData
maxRatio <- suppressWarnings(max(apply(ratioPassive[,siteColumns], 1, max, na.rm=TRUE) ))
  

## ----echo=TRUE, eval=TRUE, warning=FALSE----------------------------------------------------------

siteColumnsIndex <- grep("site", names(passiveData))
maxRatioBySite <- data.frame(site=names(passiveData)[siteColumnsIndex],
                             ratio_perc=rep(NA, length(siteColumnsIndex)),
                             chemical = rep("", length(siteColumnsIndex)),
                             endpoint =  rep("", length(siteColumnsIndex)),
                             stringsAsFactors=FALSE)

for(i in names(passiveData)[siteColumnsIndex]){ # i = site columns
  oneSite <- passiveData[,c(commonColumns, i)]
  names(oneSite)[names(oneSite) == i] <- 'value'
  casnRow <- setNames(1:nrow(oneSite),oneSite$CAS) # casnRow gets index of CAS in oneSite
  ratioPassive <- oneSite[casnRow[endPoint$casn],"value"] / endPointData
  
  maxRatio <- max(apply(ratioPassive[,siteColumnsIndex], 1, max, na.rm=TRUE) )*100
  maxRatioBySite[maxRatioBySite$site == i,"ratio_perc"] <- maxRatio
  
  if(is.finite(maxRatio)){
    maxIndexChemical <- which.max(apply(ratioPassive[,siteColumnsIndex], 1, max, na.rm=TRUE) )
    maxIndexEndpoint <- apply(ratioPassive[,siteColumnsIndex], 1, which.max)[[maxIndexChemical]]
    maxRatioBySite[maxRatioBySite$site == i,"chemical"] <- passiveData$Chemical[maxIndexChemical]
    maxRatioBySite[maxRatioBySite$site == i,"endpoint"] <- names(endPoint)[maxIndexEndpoint]
  }
}

maxRatioBySite <- arrange(maxRatioBySite, desc(ratio_perc)) %>%
  filter(ratio_perc > 1E-4) %>%
  rename("Ratio [%]"=ratio_perc, 
         "Chemical"=chemical, 
         "End Point" = endpoint) 
kable(maxRatioBySite, digits=2, row.names = FALSE)


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
packagePath <- system.file("extdata", package="toxEval")
filePath <- file.path(packagePath, "pCodeInfo.RData")
load(file=filePath)

pCodeInfo <- pCodeInfo[pCodeInfo$casrn != "", ]

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
  filter(EAR > 1)

kable(detLevels, digits = 3)


## -------------------------------------------------------------------------------------------------
qualColumns <- grep("qualifier", names(waterSamples))
valColumns <- grep("valueToUse", names(waterSamples))

waterData <- waterSamples[,valColumns]
waterData[waterSamples[,qualColumns] == "<"] <- NA



## ----echo=FALSE, warning=FALSE--------------------------------------------------------------------

dataSummary  <- data.frame(pcode=names(waterData),
                         maxValue=apply(waterData, 2, max, na.rm=TRUE),
                         row.names=NULL, stringsAsFactors=FALSE)  %>%
  mutate(pCode=sapply(strsplit(pcode, "_"), function(x) x[2])) %>%
  left_join(pCodeInfo[c("parameter_cd","casrn")], by=c("pCode"="parameter_cd")) %>%
  select(casrn, pCode, maxValue) %>%
  right_join(maxMinSummary, by=c("casrn"="casn")) %>%
  mutate(EAR=maxValue/minEndPoint) %>%
  arrange(desc(EAR)) %>%
  filter(EAR > 1)
  

kable(dataSummary, digits=3)


