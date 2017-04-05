## ----setup, include=FALSE---------------------------------
library(knitr)
library(rmarkdown)
options(continue=" ")
options(width=60)
knitr::opts_chunk$set(echo = TRUE)

## ----warning=FALSE, message=FALSE, fig.width=8, fig.height=10----
library(readxl)
library(toxEval)
library(dplyr)
library(tidyr)
library(ggplot2)

path_to_tox <-  system.file("extdata", package="toxEval")
file_name <- "OWC_data_fromSup.xlsx"

full_path <- file.path(path_to_tox, file_name)

chem_data <- read_excel(full_path, sheet = "Data")
chem_info <- read_excel(full_path, sheet = "Chemicals") 
chem_site <- read_excel(full_path, sheet = "Sites")

#Remove DEET:
chem_data <- chem_data[chem_data$CAS != "134-62-3",] 
chem_info <- chem_info[chem_info$CAS != "134-62-3",] 
#Trim names for graph:
chem_info$Class[chem_info$Class == "Human Drug, Non Prescription"] <- "Human Drug"
chem_info$Class[chem_info$Class == "Antimicrobial Disinfectant"] <- "Antimicrobial"
chem_info$Class[chem_info$Class == "Detergent Metabolites"] <- "Detergent"

ACClong <- get_ACC(chem_info$CAS)
ACClong <- remove_flags(ACClong)

cleaned_ep <- clean_endPoint_info(endPointInfo)
filtered_ep <- filter_groups(cleaned_ep)

chemicalSummary <- get_chemical_summary(ACClong,
                                        filtered_ep,
                                        chem_data, 
                                        chem_site, 
                                        chem_info)


# Make a "summary" of EAR's using Water Quality Guidelines instead of toxCast:
guideline_sum <- chem_info %>%
  gather(endPoint, WQ_value, -CAS) %>%
  filter(WQ_value != "-") %>%
  mutate(WQ_value = as.numeric(WQ_value)) %>%
  right_join(chem_data, by="CAS")

WQ <-  guideline_sum %>%
  filter(endPoint %in% c("AqT_EPA_acute","AqT_EPA_chronic",
                         "AqT_other_chronic","AqT_other_acute")) %>%
  mutate(EAR = Value/WQ_value) %>%
  group_by(SiteID,`Sample Date`,CAS) %>%
  summarise(sumEAR=sum(EAR)) %>%
  data.frame() %>%
  group_by(SiteID, CAS) %>%
  summarise(maxEAR=max(sumEAR)) %>%
  data.frame() %>%
  mutate(guide_up = "Water Quality Guideline") %>%
  mutate(guide_side = "Traditional")

EEQ <- guideline_sum %>%
  filter(endPoint %in% c("EEF_avg_in.vitro","EEF_max_in.vitro_or_in.vivo")) %>%
  mutate(EAR = Value*WQ_value*1000) %>%
  group_by(SiteID,`Sample Date`,CAS) %>%
  summarise(sumEAR=sum(EAR)) %>%
  data.frame() %>%
  group_by(SiteID, CAS) %>%
  summarise(maxEAR=max(sumEAR)) %>%
  data.frame() %>%
  mutate(guide_up = "EEQ") %>%
  mutate(guide_side = "Traditional")

# Make a "summary" of EAR's using toxCast:

toxCast <- chemicalSummary %>%
  select(-Bio_category, -shortName, -chnm) %>%
  group_by(site,date,casrn, Class) %>%
  summarise(sumEAR=sum(EAR)) %>%
  data.frame() %>%
  group_by(site, casrn, Class) %>%
  summarise(maxEAR=max(sumEAR)) %>%
  data.frame() %>%
  rename(SiteID = site,
         CAS = casrn) %>%
  mutate(guide_side = "ToxCast") %>%
  left_join(select(chem_info, CAS, `Chemical Name`), by="CAS") 

# We need extra rows for comparing tox with WQG and EEQ:
tox_WQ <- toxCast %>%
  mutate(guide_up = "Water Quality Guideline")

tox_EEQ <- filter(toxCast, CAS %in% unique(EEQ$CAS)) %>%
  mutate(guide_up = "EEQ")

cbValues <- c("#DCDA4B","#999999","#00FFFF","#CEA226","#CC79A7","#4E26CE",
              "#FFFF00","#78C15A","#79AEAE","#FF0000","#00FF00","#B1611D",
              "#FFA500","#F4426e")

order_Class <- toxCast %>%
  group_by(Class,`Chemical Name`) %>%
  summarise(median = median(maxEAR[maxEAR != 0])) %>%
  data.frame() %>%
  arrange(desc(median)) %>%
  filter(!duplicated(Class)) %>%
  arrange(median) 

EEQ <- EEQ %>%
  left_join(select(chem_info, CAS, `Chemical Name`, Class), by="CAS")
WQ <- WQ %>%
  left_join(select(chem_info, CAS, `Chemical Name`, Class), by="CAS")


graphData <- bind_rows(tox_WQ, tox_EEQ, EEQ, WQ) 

#Primary ordering needs to be the highest -> lowest class, then order by median
orderChem <- bind_rows(tox_WQ, 
                       filter(EEQ, !(CAS %in% unique(tox_WQ$CAS))),
                       filter(WQ, !(CAS %in% unique(tox_WQ$CAS)))) %>%
  group_by(`Chemical Name`,Class) %>%
  summarise(median = quantile(maxEAR[maxEAR != 0],0.5)) %>%
  data.frame() %>%
  mutate(Class = factor(Class, levels=order_Class$Class)) %>%
  arrange(Class, !is.na(median), median)

orderedLevels <- as.character(orderChem$Chemical.Name) 


graphData <-graphData %>%
  mutate(Class = factor(Class, levels=rev(order_Class$Class)),
         `Chemical Name` = factor(`Chemical Name`, levels=orderedLevels),
         guide_side = factor(guide_side),
         guide_up = factor(guide_up, levels = c("Water Quality Guideline","EEQ")))  
  
levels(graphData$guide_side) <- c("ToxCast\nMaximum EAR Per Site",
                                  "Traditional\nMaximum Quotient Per Site")

#Adding counts to the side:
countNonZero <- graphData %>%
  select(SiteID, `Chemical Name`,guide_side,guide_up, maxEAR) %>%
  group_by(SiteID, `Chemical Name`,guide_side,guide_up) %>%
  summarise(meanEAR = mean(maxEAR, na.rm=TRUE)) %>%
  group_by(`Chemical Name`,guide_side,guide_up) %>%
  summarise(nonZero = as.character(sum(meanEAR>0))) %>%
  data.frame() %>%
  select(Chemical.Name, guide_up, nonZero) %>%
  distinct() %>%
  mutate(guide_side = factor("ToxCast\nMaximum EAR Per Site", 
                             levels = levels(graphData$guide_side)),
         guide_up = factor(guide_up, levels = levels(graphData$guide_up)),
         `Chemical Name` = factor(Chemical.Name, 
                                  levels = levels(graphData$`Chemical Name`))) 

# WQ: Astricts to chemicals with no endpoints:
astrictData_WQ <- countNonZero %>%
  mutate(guide_side = factor("Traditional\nMaximum Quotient Per Site", 
                             levels = levels(graphData$guide_side))) %>%
  filter(guide_up == "Water Quality Guideline") %>%
  mutate(nonZero = "*") %>%
  filter(!(`Chemical Name` %in% unique(WQ$`Chemical Name`)))

# EEQ: Astricts to chemicals with no endpoints:
astrictData_EEQ <- countNonZero %>%
  mutate(guide_side = factor("ToxCast\nMaximum EAR Per Site", 
                             levels = levels(graphData$guide_side))) %>%
  filter(guide_up == "EEQ") %>%
  mutate(nonZero = "*") %>%
  filter(!(`Chemical Name` %in% unique(tox_EEQ$`Chemical Name`)))

# Label upper right corner for each facet (probably an easier way...):
textData <- select(graphData, guide_up, guide_side) %>%
  distinct() %>%
  mutate(textExplain = c("A","B","C","D"),
         y = c(10,10,100,100),
         `Chemical Name` = factor(rep("4-Nonylphenol",4), levels = levels(graphData$`Chemical Name`)))
  
toxPlot_All <- ggplot(data=graphData) +
  scale_y_log10(labels=fancyNumbers)  +
  geom_boxplot(aes(x=`Chemical Name`, y=maxEAR, fill=Class),
               lwd=0.1,outlier.size=1) +
  facet_grid(guide_up ~ guide_side, scales = "free", space = "free") +
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
  geom_text(data=countNonZero, aes(x= `Chemical Name`, label = nonZero, y=ymin), size=2.5) +
  geom_text(data = textData, aes(x=`Chemical Name`, label=textExplain, y=y),
            size = 3) +
  geom_text(data = astrictData_WQ, aes(x=`Chemical Name`, label=nonZero, y=10^-5),
            size=5, vjust = 0.70) +
  geom_text(data = astrictData_EEQ, aes(x=`Chemical Name`, label=nonZero, y=3.3*ymin),
            size=5, vjust = 0.70)

toxPlot_All_withLabels



