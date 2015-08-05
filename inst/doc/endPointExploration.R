## ----message=FALSE, echo=TRUE, results='asis'----------------------------
library(toxEval)
library(DT)
library(dplyr)
library(tidyr)

wData <- wData
pCodeInfo <- pCodeInfo

packagePath <- system.file("extdata", package="toxEval")
filePath <- file.path(packagePath, "stationINFO.RData")
load(file=filePath)

endPoint <- endPointToxCreate(pCodeInfo)
chemicalSummary <- chemSummBasic(wData,pCodeInfo,endPoint)
endPointInfo <- endPointInfo
# epSumm <- endPointSumm(chemicalSummary, endPointInfo=endPointInfo)

siteKey <- setNames(stationINFO$shortName, stationINFO$fullSiteID)
summary <- siteSumm(chemicalSummary,siteKey)


## ------------------------------------------------------------------------
group <- "transcription factor"
groupCol <- "intended_target_type_sub"
chemicalSummary.site <- "site"
chemicalSummary.Key <- "endPoint"
chemicalSummary.Date <- "date"
funcToSummerise <- "sumEAR"

endPointInfoSub <- select_(endPointInfo, "assay_component_endpoint_name", groupCol) %>%
  unique()

chemicalSummary <- chemicalSummary %>%
  filter(endPoint %in% endPointInfoSub$assay_component_endpoint_name ) %>%
  left_join(endPointInfoSub,
            by=c("endPoint"="assay_component_endpoint_name"))

groupSumm <- chemicalSummary %>%
  group_by(site) %>%
  summarise(nChem = length(unique(chnm)),
            nEndPoints = length(unique(endPoint))) %>%
      select(nChem,nEndPoints) %>%
      mutate(nChem = median(nChem),
             nEndPoints = median(nEndPoints)) %>%
      unique()

endpointSummary <- chemicalSummary %>%
  filter_(paste0(groupCol," == '", group, "'")) %>%
  select_("hits","EAR","endPoint",
          chemicalSummary.site,chemicalSummary.Date,groupCol) %>%
  group_by_(chemicalSummary.site, chemicalSummary.Date) %>%
  summarise(sumEAR = sum(EAR),
            nHits = sum(hits)) 

statsOfSum <- endpointSummary%>%
  group_by_(chemicalSummary.site) %>%
  summarise_(meanEAR = paste0("mean(",funcToSummerise,")"),
             # medianEAR = paste0("median(",funcToSummerise,")"),
             maxEAR = paste0("max(",funcToSummerise,")"),
             sumHits = "sum(nHits)",
             nSamples = "n()")%>%
  data.frame()%>%
  mutate(site = siteKey[site])
  




datatable(statsOfSum, rownames = FALSE, 
          caption=paste("Summary of", funcToSummerise, "of",group))  %>% 
    formatRound(c("meanEAR","maxEAR","sumHits","nSamples"), 
                digits = 1)



## ----fig.width=7, fig.height=7-------------------------------------------
library(ggplot2)

siteLimits <- stationINFO %>%
  mutate(lakeCat = factor(Lake, 
            levels=c("Lake Superior","Lake Michigan",
                     "Lake Huron", "St. Lawrence River",
                     "Detroit River and Lake St. Clair","Lake Erie","Lake Ontario"))) %>%
  arrange(lakeCat) %>%
  mutate(lakeColor = c("red","black","green","brown","brown","brown","blue")[as.numeric(lakeCat)] ) %>%
  filter(fullSiteID %in% unique(endpointSummary$site))
  

endPointSummBP <- endpointSummary %>%
  data.frame()%>%
  mutate(site = siteKey[site]) %>%
  mutate(site = factor(site, levels=siteLimits$Station.shortname)) 

sToxWS <- ggplot(endPointSummBP, aes(x=site, y=sumEAR)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.25, colour=siteLimits$lakeColor), 
        legend.position = "none")+
  scale_x_discrete(limits=siteLimits$Station.shortname) +
  # scale_y_log10(limits=c(0.03,5000)) +
  # labs(x = "Site")  +
  geom_text(data=data.frame(), aes(x=c(5, 18,31,45,56),
                                   y=-.5, label=c("Superior","Michigan","Huron","Erie","Ontario")),
                                   colour=factor(c("red","black","green","brown","blue"),
                                                 levels=c("red","black","green","brown","blue")), size=3) 

sToxWS

