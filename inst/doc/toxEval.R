## ----setup, include=FALSE-------------------------------------------------------------------------
library(xtable)
options(continue=" ")
options(width=100)
library(knitr)

## -------------------------------------------------------------------------------------------------

library("webchem")

inchk <- cts_convert(query = '3380-34-5', from = 'CAS', to = 'inchikey')
info <- cts_compinfo(inchikey = inchk)
info$molweight


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

# #Let's remove rows that now have all NA's:
# passiveData <- passiveData[rowSums(is.na(passiveData[,siteColumns]))!=length(siteColumns) ,]


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
  

## ----echo=TRUE, eval=TRUE-------------------------------------------------------------------------

packagePath <- system.file("extdata", package="toxEval")
filePath <- file.path(packagePath, "waterSamples.RData")
load(file=filePath)

head(waterSamples[,1:10])


## ----message=FALSE--------------------------------------------------------------------------------
library(dataRetrieval)

waterSamplePCodes <- names(waterSamples)[grep("valueToUse", names(waterSamples))]
waterSamplePCodes <- sapply(strsplit(waterSamplePCodes, "_"), function(x) x[2])

pCodeInfo <- readNWISpCode(waterSamplePCodes)

unique(pCodeInfo$parameter_units)


## -------------------------------------------------------------------------------------------------

AC50 <- left_join(AC50gain, passiveData[,c("CAS", "Units", "mlWt")],
                   by= c("casn" = "CAS")) %>%
  filter(!is.na(Units)) %>%
  rename(desiredUnits = Units) %>%
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


