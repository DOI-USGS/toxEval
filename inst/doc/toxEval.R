## ----setup, include=FALSE-------------------------------------------------------------------------
library(xtable)
options(continue=" ")
options(width=100)
library(knitr)

## -------------------------------------------------------------------------------------------------

library("webchem")

molweight <- cir_query('3380-34-5', "mw")
molweight


## ----echo=TRUE, eval=TRUE-------------------------------------------------------------------------

packagePath <- system.file("extdata", package="toxEval")
filePath <- file.path(packagePath, "passiveData.RData")
load(file=filePath)

head(passiveData[,1:8])

# Unique units:
unique(passiveData$Units)


## ----warning=FALSE--------------------------------------------------------------------------------
siteColumns <- grep("site",names(passiveData))
passiveData[,siteColumns] <- sapply(passiveData[,siteColumns], function(x) as.numeric(x))


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
kable(head(AC50))


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

maxMeasure <- apply(passiveData[,siteColumns], 1, max, na.rm=TRUE) 
minMeasure <- apply(passiveData[,siteColumns], 1, min, na.rm=TRUE)
dataSummary <- passiveData[,-siteColumns]
dataSummary$Max <- maxMeasure
dataSummary$Min <- minMeasure

maxMinSummary <- left_join(maxMinSummary,dataSummary, 
                           by=c("casn"="CAS", "Units"="Units"))


## -------------------------------------------------------------------------------------------------
sum(maxMinSummary$MLD > maxMinSummary$maxEndPoint)


## -------------------------------------------------------------------------------------------------
sum(maxMinSummary$Max > maxMinSummary$minEndPoint)


## ----echo=FALSE-----------------------------------------------------------------------------------
kable(head(maxMinSummary[maxMinSummary$Max > maxMinSummary$minEndPoint,c(1:5,7:8,10:11)]), digits=3)


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
  

## ----echo=TRUE, eval=FALSE------------------------------------------------------------------------
#  
#  siteColumnsIndex <- grep("site", names(passiveData))
#  maxRatioBySite <- data.frame(site=names(passiveData)[siteColumnsIndex],
#                               ratio=rep(NA, length(siteColumnsIndex)),
#                               chemical = rep("", length(siteColumnsIndex)),
#                               endpoint =  rep("", length(siteColumnsIndex))
#  
#  for(i in names(passiveData)[siteColumnsIndex]){ # i = site columns
#    oneSite <- passiveData[,c(commonColumns, i)]
#    names(oneSite)[names(oneSite) == i] <- 'value'
#    casnRow <- setNames(1:nrow(oneSite),oneSite$CAS) # casnRow gets index of CAS in oneSite
#    ratioPassive <- oneSite[casnRow[endPoint$casn],"value"] / endPointData
#  
#    maxRatio <- suppressWarnings(max(apply(ratioPassive[,siteColumnsIndex], 1, max, na.rm=TRUE) ))
#  
#  #   maxIndexChemical <- suppressWarnings(which.max(apply(ratioPassive[,siteColumnsIndex], 1, max, na.rm=TRUE) ))
#  #
#  #   maxIndexEndpoint <- apply(ratioPassive[,siteColumnsIndex], 1, which.max)[[maxIndexChemical]]
#  
#    maxRatioBySite[maxRatioBySite$site == i,"ratio"] <- maxRatio
#  #   maxRatioBySite[maxRatioBySite$site == i,"chemical"] <- passiveData$Chemical[maxIndexChemical]
#  #   maxRatioBySite[maxRatioBySite$site == i,"endpoint"] <- names(endPoint)[maxIndexEndpoint]
#  
#  }
#  
#  #So max ratio is:
#  max(maxRatioBySite$ratio)
#  
#  # #At site:
#  # maxRatioBySite$site[maxRatioBySite$ratio == max(maxRatioBySite$ratio)]
#  #
#  # #The chemical that gave the max ratio:
#  
#  

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

# Check if any changes over time:
all(apply(detectionLimits, 2, function(x) length(unique(x))==1))

#It does!


