## ----setup, include=FALSE-------------------------------------------------------------------------
options(continue=" ")
options(width=100)
library(knitr)
library(rmarkdown)
library(dplyr)
library(reshape2)
library(DT)
library(ggplot2)


## ----echo=FALSE, eval=TRUE, message=FALSE---------------------------------------------------------
library(toxEval)

packagePath <- system.file("extdata", package="toxEval")
#Site Info:
filePath <- file.path(packagePath, "stationINFO.RData")
load(file=filePath)
#parameter info:
filePath <- file.path(packagePath, "pCodeInfo.RData")
load(file=filePath)


pCodeInfo <- pCodeInfo[pCodeInfo$casrn != "", ]



## ----message=FALSE, echo=FALSE--------------------------------------------------------------------

# AC50 data provided in the toxEval package:
AC50gain <- select(pCodeInfo, srsname, casrn, AqT_EPA_acute, 
                   AqT_EPA_chronic, AqT_other_acute, AqT_other_chronic)


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
#  pCodeInfo$mlWt[pCodeInfo$casrn != ""] <- sapply(pCodeInfo$casrn[
#    pCodeInfo$casrn != ""],
#                  function(x) cir_query(x, "mw", first = TRUE))
#  pCodeInfo$mlWt <- as.numeric(pCodeInfo$mlWt)
#  

## ----message=FALSE, echo=FALSE--------------------------------------------------------------------

waterSamplePCodes <- names(waterSamples)[grep("valueToUse", names(waterSamples))]
waterSamplePCodes <- sapply(strsplit(waterSamplePCodes, "_"), function(x) x[2])




## -------------------------------------------------------------------------------------------------

AC50 <- left_join(AC50gain, pCodeInfo[pCodeInfo$parameter_cd %in% waterSamplePCodes,c("casrn", "parameter_units", "mlWt")],
                   by= c("casrn"="casrn")) %>%
  filter(!is.na(parameter_units)) %>%
  rename(desiredUnits = parameter_units) %>%
  mutate(conversion = 1) %>%
  select(casrn, srsname, desiredUnits, mlWt, conversion) 


## -------------------------------------------------------------------------------------------------
AC50Converted <- left_join(AC50, AC50gain, by = c("casrn", "srsname")) %>%
  rename(casn=casrn, chnm=srsname)

infoColumns <- c("casn", "chnm", "desiredUnits","mlWt", "conversion")

endPointData <- AC50Converted[,!(names(AC50Converted) %in% infoColumns)]

endPointData <- endPointData * AC50Converted$conversion

AC50 <- rename(AC50, casn=casrn, chnm=srsname)

endPoint <- cbind(AC50, data.frame(endPointData))
endPoint <- rename(endPoint, Units=desiredUnits)
# endPoint <- na.omit(endPoint)

## ----warning=FALSE--------------------------------------------------------------------------------
maxEndPoints <- apply(endPointData, 1, max, na.rm=TRUE) 
minEndPoints <- apply(endPointData, 1, min, na.rm=TRUE)

maxMinSummary <- cbind(endPoint[,c("chnm", "casn", "Units")], 
                       maxEndPoint=maxEndPoints, 
                       minEndPoint=minEndPoints)


## -------------------------------------------------------------------------------------------------
pCodeSummary <- select(pCodeInfo, srsname, casrn, parameter_units, mlWt) %>%
  filter(pCodeInfo$parameter_cd %in% waterSamplePCodes)

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
  left_join(pCodeInfo[c("parameter_cd","casrn","class")], 
            by=c("pCode"="parameter_cd")) %>%
  select(casrn, class, minDetLevel) %>%
  right_join(maxMinSummary, by=c("casrn"="casn")) %>%
  mutate(EAR=minDetLevel/minEndPoint) %>%
  arrange(desc(EAR)) %>%
  filter(EAR > 0.1) %>%
  rename(Chemical=chnm) %>%
  select(Chemical, EAR, minDetLevel, minEndPoint, class, Units)

datatable(detLevels, rownames = FALSE,options = list(pageLength = 25))  %>% 
    formatRound(c("EAR","minDetLevel", "minEndPoint"), digits = 3)


## -------------------------------------------------------------------------------------------------
qualColumns <- grep("qualifier", names(waterSamples))
valColumns <- grep("valueToUse", names(waterSamples))

waterData <- waterSamples[,valColumns]
waterData[waterSamples[,qualColumns] == "<"] <- 0


## ----message=FALSE, results='asis', echo=FALSE----------------------------------------------------
newSiteKey <- setNames(stationINFO$shortName, stationINFO$fullSiteID)

wData <- cbind(waterSamples[,1:2],waterData)
valColumns <- grep("valueToUse",names(wData))

chemicalSummary <- wData %>%
  melt(id.vars=c("ActivityStartDateGiven","site")) %>%
  mutate(variable=as.character(variable)) %>%
  filter(!is.na(value)) %>%
  rename(measuredValue=value, pCode=variable, date=ActivityStartDateGiven) %>%
  mutate(pCode=sapply(strsplit(pCode, "_"), function(x) x[2])) %>%
  left_join(pCodeInfo[c("parameter_cd","casrn","class","srsname")], by=c("pCode"="parameter_cd")) %>%
  select(srsname,casrn, class, measuredValue, site, date) %>%
  right_join(endPoint, by=c("casrn"="casn")) %>%
  select(-mlWt, -conversion, -casrn,  -Units, -srsname) %>%
  melt(id.vars = c("class", "site", "measuredValue","chnm","date")) %>% 
  filter(!is.na(value)) %>%
  mutate(variable=as.character(variable)) %>%
  rename(endPointValue=value, endPoint=variable) %>%
  mutate(EAR=measuredValue/endPointValue) 


chemSum1 <- chemicalSummary %>%
  group_by(chnm,  class, date, site) %>%
  summarize(hits= as.numeric(any(EAR > 0.1))) %>%
  group_by(chnm,  class, site) %>%
  summarize(hits=as.numeric(any(hits > 0)),
            counts=n(),
            freq=hits/counts) %>%
  group_by(chnm, class) %>%
  summarize(freq=sum(hits)/n()) 

chemSum2 <- chemicalSummary %>%
  group_by(chnm,  class) %>%  
  summarize(maxEAR=max(EAR), 
            nSamples=nrow(unique(data.frame(date,site))),
            nEndPoints=length(unique(endPoint[EAR > 0.1]))) %>%
  left_join(chemSum1, by=c("chnm","class")) %>%
  data.frame %>%
  arrange(desc(freq)) %>%
  select(chnm, class, freq, maxEAR, nEndPoints, nSamples)


  datatable(chemSum2, rownames = FALSE, 
            options = list(pageLength = 10)) %>% 
    formatRound(c("maxEAR", "freq"), digits = 4)



## ----message=FALSE, results='asis', echo=FALSE----------------------------------------------------

chemSum2 <- chemicalSummary  %>%
  group_by(chnm, endPoint, class) %>% 
  summarize(maxEAR=max(EAR), 
#there are some sites that have sample times at exactly the same time!:
            nSamples=nrow(unique(data.frame(date,site))), 
            freq=sum(EAR > 0.1)/nSamples, 
            nSites=length(unique(site[EAR>0.1]))) %>% 
  select(chnm, maxEAR, freq, nSamples, nSites, endPoint, class) %>%
  data.frame %>%
  filter(maxEAR > 0.1) %>%
  arrange(desc(maxEAR))
  

  datatable(unique(chemSum2), rownames = FALSE, 
            options = list(pageLength = 10)) %>% 
    formatRound(c("maxEAR", "freq"), digits = 2)


## ----message=FALSE, results='asis', echo=FALSE----------------------------------------------------

siteSummary <- wData %>%
  melt(id.vars=c("ActivityStartDateGiven","site")) %>%
  mutate(variable=as.character(variable)) %>%
  filter(!is.na(value)) %>%
  rename(measuredValue=value, pCode=variable, date=ActivityStartDateGiven) %>%
  mutate(pCode=sapply(strsplit(pCode, "_"), function(x) x[2])) %>%
  left_join(pCodeInfo[c("parameter_cd","casrn","class","srsname")], by=c("pCode"="parameter_cd")) %>%
  select(srsname,casrn, class, measuredValue, site, date) %>%
  right_join(endPoint, by=c("casrn"="casn")) %>%
  select(-mlWt, -conversion, -casrn,  -Units, -srsname) %>%
  melt(id.vars = c("class", "site", "date","measuredValue","chnm")) %>% 
  filter(!is.na(value)) %>%
  mutate(variable=as.character(variable)) %>%
  rename(endPointValue=value, endPoint=variable) %>%
  mutate(EAR=measuredValue/endPointValue) %>%
  mutate(site=newSiteKey[site]) %>%
  select(site, chnm, EAR, endPoint, class, date) %>% 
  arrange(site, chnm, EAR)

summ1 <- siteSummary %>%
  group_by(site, date) %>%
  summarize(hits= as.numeric(any(EAR > 0.1))) %>%
  group_by(site) %>%
  summarize(freq=sum(hits)/length(unique(date)),
            nSamples=length(unique(date))) 

summ2 <-  siteSummary %>%
  group_by(site) %>%
  summarize(nChem = length(unique(chnm[EAR > 0.1])),
            nEndPoints = length(unique(endPoint[EAR > 0.1])),
            maxEAR = max(EAR))


summary <- left_join(summ1, summ2, by="site") %>%
  arrange(desc(maxEAR))

datatable(summary, rownames = FALSE) %>% 
#           options = list(pageLength = 50)) %>% 
  formatRound(c("freq","maxEAR"), digits = 2)


## ----message=FALSE, results='asis', echo=FALSE----------------------------------------------------

for(i in summary$site[!is.na(summary$site)]){
  
  oneSite <- wData %>%
    mutate(site = newSiteKey[site]) %>%
    filter(site == i) 
  
  oneSiteLong <- select(oneSite, -site) %>%
    melt(id.vars=c("ActivityStartDateGiven")) %>%
    mutate(variable=as.character(variable)) %>%
    filter(!is.na(value)) %>%
    rename(measuredValue=value, pCode=variable, date=ActivityStartDateGiven) %>%
    mutate(pCode=sapply(strsplit(pCode, "_"), function(x) x[2])) %>%
    left_join(pCodeInfo[c("parameter_cd","casrn","class","srsname")], by=c("pCode"="parameter_cd")) %>%
    select(srsname,casrn, class, measuredValue, date) %>%
    right_join(endPoint, by=c("casrn"="casn")) %>%
    select(-mlWt, -conversion, -casrn,  -Units, -srsname) %>%
    melt(id.vars = c("class", "date","measuredValue","chnm")) %>% 
    filter(!is.na(value)) %>%
    mutate(variable=as.character(variable)) %>%
    rename(endPointValue=value, endPoint=variable) %>%
    mutate(EAR=measuredValue/endPointValue) %>%
  #   filter(EAR > 0.1) %>%
    select(chnm, EAR, endPoint, class, date) %>% 
    arrange(chnm, EAR) %>%
    group_by(chnm, endPoint, class) %>%
    summarize(hits=sum(EAR > 0.10),
              nSamples=length(unique(date)),
              freq=hits/nSamples) %>%
    data.frame()%>%
    filter(freq>0) %>%
    arrange(chnm, desc(hits))  

  if(nrow(oneSiteLong) > 0){
    cat("\n\n###", i, "\n")
    print(kable(oneSiteLong, digits=3,caption = i, row.names = FALSE))
  }
  

}


## ----fig.width=10, fig.height=10, warning=FALSE, echo=FALSE, fig.cap="Colors represent different endpoints, dots represent measured water sample values"----

# maxY <- 10
# data <- cbind(waterSamples[,1:2],waterData) %>%
#   melt(id.vars=c("ActivityStartDateGiven","site")) %>%
#   mutate(variable=as.character(variable)) %>%
#   filter(!is.na(value)) %>%
#   rename(measuredValue=value, pCode=variable) %>%
#   mutate(pCode=sapply(strsplit(pCode, "_"), function(x) x[2])) %>%
#   left_join(pCodeInfo[c("parameter_cd","casrn","class","srsname")], by=c("pCode"="parameter_cd")) %>%
#   select(srsname,casrn, class, measuredValue, site) %>%
#   right_join(endPoint, by=c("casrn"="casn")) %>%
#   select(-mlWt, -conversion, -casrn,  -Units, -srsname) %>%
#   melt(id.vars = c("class", "site", "measuredValue","chnm")) %>% 
#   filter(!is.na(value)) %>%
#   mutate(variable=as.character(variable)) %>%
#   rename(endPointValue=value, endPoint=variable) %>%
#   mutate(endPointValue=ifelse(endPointValue <= maxY, maxY, endPointValue)) 
# 
# # ggplot(ep[1:50,], aes(variable)) + geom_bar() + facet_wrap(~ chnm)
# # ggplot(ep[1:50,], aes(variable, fill=chnm)) + geom_bar() # + scale_fill_brewer()
# 
# ggplot(data, aes(x = chnm, y = endPointValue, fill=endPoint)) +
#     geom_bar(stat='identity',guide=FALSE,colour="grey") + 
#     ylim(0,10) +   
#     theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.25), legend.position = "none") +
#     labs(x = "Chemical", y="Concentration") + 
#     geom_point(aes(x=chnm, y=measuredValue))




## ----fig.width=10, fig.height=7, warning=FALSE, echo=FALSE, fig.cap="Boxplot of all EARS over 0.1"----

# data2 <- cbind(waterSamples[,1:2],waterData) %>%
#   melt(id.vars=c("ActivityStartDateGiven","site")) %>%
#   mutate(variable=as.character(variable)) %>%
#   filter(!is.na(value)) %>%
#   rename(measuredValue=value, pCode=variable) %>%
#   mutate(pCode=sapply(strsplit(pCode, "_"), function(x) x[2])) %>%
#   left_join(pCodeInfo[c("parameter_cd","casrn")], by=c("pCode"="parameter_cd")) %>%
#   select(casrn, measuredValue, site) %>%
#   right_join(endPoint, by=c("casrn"="casn")) %>%
#   select(-mlWt, -conversion, -casrn,  -Units) %>%
#   melt(id.vars = c("site", "measuredValue","chnm")) %>% 
#   filter(!is.na(value)) %>%
#   mutate(variable=as.character(variable)) %>%
#   rename(endPointValue=value, endPoint=variable) %>%
#   mutate(EAR = measuredValue/endPointValue) %>%
#   filter(EAR >= 0.1)
# 
# ggplot(data2, aes(x = chnm, y = EAR)) +
#     geom_boxplot() + 
#     ylim(0,3) +    
#     theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.25), legend.position = "none") +
#     labs(x = "Chemical", y="EAR")



## ----fig.width=10, fig.height=7, warning=FALSE, echo=FALSE, fig.cap="Boxplot of all max EARS per site per chemical (filtered over 0.1)"----

# data3 <- group_by(data2, site, chnm) %>%
#   mutate(maxEAR = max(EAR))
# 
# ggplot(data3, aes(x = chnm, y = maxEAR)) +
#     geom_boxplot() + 
#     ylim(0,3) +    
#     theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.25), legend.position = "none") +
#     labs(x = "Chemical", y="maxEAR")+ 
#     geom_point(aes(x=chnm, y=EAR, color="red"))


## ----fig.width=10, fig.height=7, warning=FALSE, echo=FALSE, fig.cap="Boxplot of  EARS per site (filtered over 0.1, max per endpoint)"----

# data4 <- mutate(data2, site = newSiteKey[site]) %>%
#   group_by(site, chnm, endPoint) %>%
#   mutate(maxEAR = max(EAR))
#   
# 
# ggplot(data4, aes(x = site, y = maxEAR)) +
#     geom_boxplot() + 
#     ylim(0,3) +    
#     theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.25), legend.position = "none") +
#     labs(x = "Site", y="maxEAR")+ 
#     geom_point(aes(x=site, y=EAR, color="red") )


