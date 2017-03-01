library(toxEval)
library(dplyr)
library(data.table)
library(tidyr)
library(ggplot2)
library(grid)
library(gridExtra)

source("getDataReady.R")
# source("getCenDataReady.R")
###########################################################################
# WQ boxplot
filePath <- file.path(pathToApp, "waterSamples.RData")
load(file=filePath)
waterSamples$site["USGS-04157005" == waterSamples$site] <- "USGS-04157000"
valColumns <- grep("valueToUse", names(waterSamples))
qualColumns <- grep("qualifier", names(waterSamples))
waterData <- waterSamples[,valColumns]
waterData[waterSamples[,qualColumns] == "<"] <- 0
wData <- cbind(waterSamples[,1:2],waterData)

wDataLong <- gather(wData, pCode, measuredValue, -ActivityStartDateGiven, -site) %>%
  rename(date = ActivityStartDateGiven) %>%
  filter(!is.na(measuredValue)) %>%
  mutate(pCode = gsub("valueToUse_", replacement = "", pCode)) %>%
  left_join(select(pCodeInfo, parameter_cd, casrn, class, chnm,
                   # EEF_avg_in.vitro,EEF_max_in.vitro_or_in.vivo,
                   AqT_EPA_acute,AqT_EPA_chronic,
                   AqT_other_acute,AqT_other_chronic),
            by=c("pCode" = "parameter_cd")) %>%
  gather(endPoint, value, 
         # EEF_avg_in.vitro,EEF_max_in.vitro_or_in.vivo,
         AqT_EPA_acute,AqT_EPA_chronic,
         AqT_other_acute,AqT_other_chronic) %>%
  mutate(EAR =  measuredValue/value) %>%
  filter(!is.na(EAR)) 

waterSamplePCodes <- unique(wDataLong$pCode)

toxCastChems <- gather(ACC, endPoint, ACC, -casn, -chnm, -flags) %>%
  filter(!is.na(ACC)) %>%
  left_join(pCodeInfo[pCodeInfo$parameter_cd %in% waterSamplePCodes,c("casrn", "parameter_units", "mlWt")],
            by= c("casn"="casrn")) %>%
  filter(!is.na(parameter_units)) %>%
  select(casn) %>%
  distinct()

wDataLong <- wDataLong %>%
  filter(casrn %in% toxCastChems$casn)

graphData.wq <- wDataLong %>%
  group_by(site,date,chnm,class) %>%
  summarise(sumEAR=sum(EAR)) %>%
  data.frame() %>%
  group_by(site, chnm,class) %>%
  summarise(meanEAR=max(sumEAR)) %>%
  data.frame() 

##################################################
# Regular toxEval stuff:
graphData <- graphData %>%
  mutate(category = as.character(category),
         guideline = "ToxCast")

graphData.full_WQ <- mutate(graphData.wq, class=as.character(class)) %>%
  rename(category = chnm) %>%
  mutate(guideline = "Traditional")

graphData.full_WQ$class[graphData.full_WQ$class == "Detergent Metabolites"] <- "Detergent"

wDataLong_EQ <- gather(wData, pCode, measuredValue, -ActivityStartDateGiven, -site) %>%
  rename(date = ActivityStartDateGiven) %>%
  filter(!is.na(measuredValue)) %>%
  mutate(pCode = gsub("valueToUse_", replacement = "", pCode)) %>%
  left_join(select(pCodeInfo, parameter_cd, casrn, class, chnm,
                   EEF_max_in.vitro_or_in.vivo),
            by=c("pCode" = "parameter_cd")) %>%
  gather(endPoint, value, EEF_max_in.vitro_or_in.vivo) %>%
  mutate(EAR =  measuredValue*value*1000) %>% # we were dividing by 0.7...not sure why
  filter(!is.na(EAR)) 

waterSamplePCodes_EQ <- unique(wDataLong_EQ$pCode)

toxCastChems_EQ <- gather(ACC, endPoint, ACC, -casn, -chnm, -flags) %>%
  filter(!is.na(ACC)) %>%
  left_join(pCodeInfo[pCodeInfo$parameter_cd %in% waterSamplePCodes_EQ,c("casrn", "parameter_units", "mlWt")],
            by= c("casn"="casrn")) %>%
  filter(!is.na(parameter_units)) %>%
  select(casn) %>%
  distinct()

# toxCastChems_EQ$chnm[toxCastChems_EQ$chnm == "4-(1,1,3,3-Tetramethylbutyl)phenol"] <- "4-tert-Octylphenol"

# wDataLong_EQ <- wDataLong_EQ %>%
#   filter(casrn %in% toxCastChems_EQ$casn)

graphData.eq <- wDataLong_EQ %>%
  group_by(site,date,chnm,class) %>%
  summarise(sumEAR=sum(EAR)) %>%
  data.frame() %>%
  group_by(site, chnm,class) %>%
  summarise(meanEAR=max(sumEAR)) %>%
  data.frame() %>%
  mutate(guideline = "Traditional") %>%
  rename(category = chnm)

subTox <- filter(graphData, category %in% graphData.eq$category) %>%
  mutate(otherThing = "EEQ")
  
# subTox$category[subTox$category == "4-(1,1,3,3-Tetramethylbutyl)phenol"] <- "4-tert-Octylphenol"

EQ <- graphData.eq %>%
  mutate(otherThing = "EEQ")

# EQ$category[EQ$category == "4-(1,1,3,3-Tetramethylbutyl)phenol"] <- "4-tert-Octylphenol"

subToxWQ <- graphData %>%
  mutate(otherThing = "Water Quality Guidelines")

WQ <- graphData.wq %>%
  rename(category = chnm) %>%
  mutate(otherThing = "Water Quality Guidelines") %>%
  mutate(guideline = "Traditional")

fullFULL <- bind_rows(subTox, EQ, subToxWQ, WQ) 
fullFULL$class[fullFULL$class == "Detergent Metabolites"] <- "Detergent"

fullFULL <- fullFULL %>%  mutate(class = factor(class, levels=rev(as.character(orderClass$class))))

fullData <- bind_rows(graphData.full_WQ, graphData) 

fullData$class[fullData$class == "Detergent Metabolites"] <- "Detergent"

orderChem <- graphData %>%#fullData %>% #not fullFULL...or just graphData....needs just tox and WQ
  group_by(category,class) %>%
  summarise(median = quantile(meanEAR[meanEAR != 0],0.5)) %>%
  data.frame() %>%
  mutate(class = factor(class, levels=orderClass$class)) %>%
  arrange(class, median)

orderedLevels <- as.character(orderChem$category)
orderedLevels <- orderedLevels[!is.na(orderedLevels)]
orderedLevels <- c(orderedLevels[1:2], "Isopropylbenzene (cumene)", 
                   orderedLevels[3:4],"Tribromomethane (bromoform)",
                   orderedLevels[5:length(orderedLevels)])
orderedLevels <- c(orderedLevels[1:46], "4-Nonylphenol diethoxylate",             
                   "4-Nonylphenol monoethoxylate",            
                   "4-tert-Octylphenol monoethoxylate (OP1EO)",
                   "4-tert-Octylphenol diethoxylate (OP2EO)", orderedLevels[47:48] )

fullFULL <- fullFULL %>%
  mutate(guideline = factor(as.character(guideline), levels=c("ToxCast","Traditional")),
         otherThing = factor(as.character(otherThing), levels = c("Water Quality Guidelines","EEQ"))) %>%
  mutate(class = factor(class, levels=rev(orderClass$class))) %>%
  mutate(category = factor(category, levels=orderedLevels)) 

fullFULL$class[fullFULL$category == "4-Nonylphenol"] <- "Detergent"

cbValues <- c("#DCDA4B","#999999","#00FFFF","#CEA226","#CC79A7","#4E26CE",
              "#FFFF00","#78C15A","#79AEAE","#FF0000","#00FF00","#B1611D",
              "#FFA500","#F4426e")

countNonZero <- fullFULL %>%
  select(site, category,guideline,otherThing, meanEAR) %>%
  group_by(site, category,guideline,otherThing) %>%
  summarise(meanEAR = mean(meanEAR, na.rm=TRUE)) %>%
  group_by(category,guideline,otherThing) %>%
  summarise(nonZero = as.character(sum(meanEAR>0))) %>%
  data.frame() %>%
  select(category, otherThing, nonZero) %>%
  distinct() %>%
  mutate(guideline = factor(c("ToxCast"), levels = levels(fullFULL$guideline)),
         otherThing = factor(otherThing, levels = levels(fullFULL$otherThing)),
         category = factor(category, levels = levels(fullFULL$category))) 

countNonZero <- countNonZero[!duplicated(countNonZero[,1:2]),]

astrictData <- countNonZero %>%
  mutate(guideline = factor(c("Traditional"), levels = levels(fullFULL$guideline))) %>%
  filter(otherThing == "Water Quality Guidelines") %>%
  mutate(nonZero = "*") %>%
  filter(!(category %in% unique(WQ$category)))

toxAst <- data.frame(category = c("4-Nonylphenol diethoxylate",             
                                  "4-Nonylphenol monoethoxylate",            
                                  "4-tert-Octylphenol monoethoxylate (OP1EO)",
                                  "4-tert-Octylphenol diethoxylate (OP2EO)",
                                  "Isopropylbenzene (cumene)",
                                  "Tribromomethane (bromoform)"),
                     otherThing = c(rep("EEQ",4), rep("Water Quality Guidelines",2))) %>%
  mutate(nonZero = "*") %>%
  mutate(category = factor(category, levels = levels(fullFULL$category))) %>%
  mutate(guideline = factor(c("ToxCast"), levels = levels(fullFULL$guideline))) 


levels(fullFULL$guideline) <- c("ToxCast\nMaximum EAR Per Site", 
                                "Traditional\nMaximum Quotient Per Site")
levels(countNonZero$guideline) <- c("ToxCast\nMaximum EAR Per Site", 
                                "Traditional\nMaximum Quotient Per Site")
levels(astrictData$guideline) <- c("ToxCast\nMaximum EAR Per Site", 
                                    "Traditional\nMaximum Quotient Per Site")
levels(toxAst$guideline) <- c("ToxCast\nMaximum EAR Per Site", 
                                   "Traditional\nMaximum Quotient Per Site")

textData <- data.frame(guideline = factor(c(rep("ToxCast\nMaximum EAR Per Site", 2),rep("Traditional\nMaximum Quotient Per Site", 2)), levels = levels(fullFULL$guideline)),
                       otherThing = factor(c("Water Quality Guidelines","EEQ",
                                             "Water Quality Guidelines","EEQ"), levels = levels(fullFULL$otherThing)),
                       category = factor(rep("4-Nonylphenol",4), levels = levels(fullFULL$category)),
                       textExplain = c("A","B","C","D"),
                       y = c(10,10,100,100))

toxPlot_All <- ggplot(data=fullFULL) +
  scale_y_log10(labels=fancyNumbers)  +
  geom_boxplot(aes(x=category, y=meanEAR, fill=class),
               lwd=0.1,outlier.size=1) +
  facet_grid(otherThing ~ guideline, scales = "free", space = "free") +
  theme_bw() +
  scale_x_discrete(drop=TRUE) +
  coord_flip() +
  theme(axis.text = element_text( color = "black"),
        axis.text.y = element_text(size=7),
        axis.title=element_blank(),
        panel.background = element_blank(),
        plot.background = element_rect(fill = "transparent",colour = NA),
        strip.background = element_rect(fill = "transparent",colour = NA),
        strip.text.y = element_blank()) +
  guides(fill=guide_legend(ncol=6)) +
  theme(legend.position="bottom",
        legend.justification = "left",
        legend.background = element_rect(fill = "transparent", colour = "transparent"),
        legend.title=element_blank(),
        legend.text = element_text(size=8),
        legend.key.height = unit(1,"line")) +
  scale_fill_manual(values = cbValues, drop=FALSE) 

ymin <- 10^-6
ymax <- ggplot_build(toxPlot_All)$layout$panel_ranges[[1]]$y.range[2]

toxPlot_All_withLabels <- toxPlot_All +
  geom_text(data=countNonZero, aes(x=category,label=nonZero, y=ymin), size=2.5) +
  geom_text(data = textData, aes(x=category, label=textExplain, y=y), 
            size = 3) +
  geom_text(data = astrictData, aes(x=category, label=nonZero, y=10^-5), 
            size=5, vjust = 0.70) +
  geom_text(data = toxAst, aes(x=category, label=nonZero, y=3.3*ymin), 
            size=5, vjust = 0.70)

toxPlot_All_withLabels

ggsave(toxPlot_All_withLabels, bg = "transparent",
       filename = "allPanels_betterEEQ.png", 
       height = 10, width = 9)
