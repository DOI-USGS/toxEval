## ----setup, include=FALSE---------------------------------
library(knitr)
library(rmarkdown)
options(continue=" ")
options(width=60)
knitr::opts_chunk$set(echo = TRUE)

## ----warning=FALSE, message=FALSE, fig.width=6, fig.height=6----
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

orderChem <- tox_WQ %>%
  group_by(`Chemical Name`,Class) %>%
  summarise(median = quantile(maxEAR[maxEAR != 0],0.5)) %>%
  data.frame() %>%
  mutate(Class = factor(Class, levels=order_Class$Class)) %>%
  arrange(Class, median)

orderedLevels <- as.character(orderChem$Chemical.Name) 
orderedLevels <- c(orderedLevels[1:2], "Isopropylbenzene (cumene)", 
                   orderedLevels[3:4],"Tribromomethane (bromoform)",
                   orderedLevels[5:length(orderedLevels)])
orderedLevels <- c(orderedLevels[1:46], "4-Nonylphenol diethoxylate",             
                   "4-Nonylphenol monoethoxylate",            
                   "4-tert-Octylphenol monoethoxylate (OP1EO)",
                   "4-tert-Octylphenol diethoxylate (OP2EO)", orderedLevels[47:48] )


graphData <-graphData %>%
  mutate(Class = factor(Class, levels=rev(order_Class$Class)),
         `Chemical Name` = factor(`Chemical Name`, levels=orderedLevels),
         guide_side = factor(guide_side),
         guide_up = factor(guide_up, levels = c("Water Quality Guideline","EEQ")))  
  
levels(graphData$guide_side) <- c("ToxCast\nMaximum EAR Per Site",
                                  "Traditional\nMaximum Quotient Per Site")

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

toxPlot_All


