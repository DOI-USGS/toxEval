source("getDataReady.R")

library(ggplot2)

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

ACC_graphData <- ACClong %>%
  filter(endPoint %in% ep$endPoint) %>%
  mutate(valueToPlot = as.numeric(value)) %>%
  left_join(select(pCodeInfo, casn=casrn, chnm, class), by="casn") %>%
  filter(chnm != "DEET") %>%
  mutate(chnm_factor = factor(chnm, levels = orderedLevels)) 

flagsShort <- c("Borderline",  "OnlyHighest",
                "GainAC50", "Biochemical")
for(i in flagsShort){
  take.out.flags <- flagDF[!flagDF[[i]],c("casn","endPoint")]
  
  ACC_graphData <- right_join(ACC_graphData, take.out.flags, 
                                by=c("casn", "endPoint"="endPoint")) %>%
    filter(!is.na(chnm))
}

ACC_graphData$class[ACC_graphData$class == "Human Drug, Non Prescription"] <- "Human Drug"
ACC_graphData$class[ACC_graphData$class == "Antimicrobial Disinfectant"] <- "Antimicrobial"
ACC_graphData$class[ACC_graphData$class == "Detergent Metabolites"] <- "Detergent"

ACC_graphData <- ACC_graphData %>%
  mutate(class = factor(class, levels=rev(orderClass$class))) 


cbValues <- c("#DCDA4B","#999999","#00FFFF","#CEA226","#CC79A7","#4E26CE",
              "#FFFF00","#78C15A","#79AEAE","#FF0000","#00FF00","#B1611D",
              "#FFA500","#F4426e")
  
ACC_graphData <- filter(ACC_graphData, chnm %in% levels(ACC_graphData$chnm_factor))

countNonZero <- ACC_graphData %>%
  select(chnm_factor, class, endPoint) %>%
  distinct() %>%
  group_by(chnm_factor, class) %>%
  summarize(nonZero = n())

ACC_plot <- ggplot(ACC_graphData) +
  scale_y_log10("ACC", labels=fancyNumbers)+
  geom_boxplot(aes(x=chnm_factor, y=valueToPlot),lwd=0.1,outlier.shape=NA) + #, fill = class,outlier.size=1
  geom_jitter(aes(x=chnm_factor, y=valueToPlot, color = class), position=position_jitter(width=.1, height=0)) +
  coord_flip() +
  theme_bw() +
  scale_x_discrete(drop=TRUE) +
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
  scale_color_manual(values = cbValues, drop=FALSE) 

ACC_plot

ymin <- 10^ggplot_build(ACC_plot)$layout$panel_ranges[[1]]$x.range[1]

ACC_plot_withLabels <- ACC_plot +
  geom_text(data=countNonZero, aes(x=chnm_factor,label=nonZero, y=ymin), size=2.5) 

ACC_plot_withLabels

ggsave(ACC_plot_withLabels, bg = "transparent",
       filename = "ACCs_withJitter.png",
       height = 7, width = 7)
