library(toxEval)
library(dplyr)
library(data.table)
library(tidyr)
library(ggplot2)
library(grid)
library(gridExtra)

source("getCenDataReady.R")

##################################################
# Regular toxEval stuff:
graphData <- chemicalSummary %>%
  select(EAR, chnm, class)

cbValues <- c("#DCDA4B","#999999","#00FFFF","#CEA226","#CC79A7","#4E26CE",
              "#FFFF00","#78C15A","#79AEAE","#FF0000","#00FF00","#B1611D",
              "#FFA500","#F4426e")

orderClass <- graphData %>%
  group_by(class,chnm) %>%
  summarise(median = median(EAR[EAR != 0])) %>%
  data.frame() %>%
  arrange(desc(median)) %>%
  filter(!duplicated(class)) %>%
  arrange(median)

orderChem <- graphData %>%
  group_by(chnm,class) %>%
  summarise(median = quantile(EAR[EAR != 0],0.5)) %>%
  data.frame() %>%
  mutate(class = factor(class, levels=orderClass$class)) %>%
  arrange(class, median)
orderedLevels <- as.character(orderChem$chnm)

graphData <- graphData %>%
  mutate(chnm = factor(chnm, levels=orderedLevels)) %>%
  mutate(class = factor(class, levels = rev(levels(orderChem$class))))

countNonZero <- chemicalSummary %>%
  select(chnm, class, endPoint) %>%
  distinct() %>%
  group_by(chnm, class) %>%
  summarize(nEndpoints = n())


toxPlot_All <- ggplot(data=graphData) +
  scale_y_log10(labels=fancyNumbers)  +
  geom_boxplot(aes(x=chnm, y=EAR, fill=class),
               lwd=0.1,outlier.size=1)  +
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

ymin <- 0.3*10^-6
ymax <- ggplot_build(toxPlot_All)$layout$panel_ranges[[1]]$y.range[2]

toxPlot_All_withLabels <- toxPlot_All +
  geom_text(data=countNonZero, aes(x=chnm,label=nEndpoints, y=ymin), size=2.5) 

# toxPlot_All_withLabels

ggsave(toxPlot_All_withLabels, bg = "transparent",
       filename = "DetectionLevels.png",
       height = 7, width = 7)

# toxPlot_All_withLabels <- toxPlot_All_withLabels + 
#   annotation_custom(
#     grob = textGrob(label = "# Endpoints", 
#                     gp = gpar(cex = 0.4)),
#     ymin = 10^-8,      # Vertical position of the textGrob
#     ymax = 10^-8,
#     xmin = 46,  # Note: The grobs are positioned outside the plot area
#     xmax = 46)
# 
# toxPlot_All_withLabels <- ggplot_gtable(ggplot_build(toxPlot_All_withLabels))
# toxPlot_All_withLabels$layout$clip[toxPlot_All_withLabels$layout$name == "panel"] <- "off"
# 
# ggsave(toxPlot_All_withLabels, #bg = "transparent",
#        filename = "DetectionLevels2.png", 
#        height = 7, width = 7)
