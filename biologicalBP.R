library(toxEval)
library(dplyr)
library(data.table)
library(tidyr)
library(ggplot2)
library(grid)
library(gridExtra)

source("getDataReady.R")

graphData <- chemicalSummary %>%
  mutate(category = choices) %>%
  group_by(site,date,category) %>%
  summarise(sumEAR=sum(EAR)) %>%
  data.frame() %>%
  group_by(site, category) %>%
  summarise(meanEAR=max(sumEAR)) %>%
  data.frame() 

graphData$class[graphData$class == "Human Drug, Non Prescription"] <- "Human Drug"
graphData$class[graphData$class == "Antimicrobial Disinfectant"] <- "Antimicrobial"
graphData$class[graphData$class == "Detergent Metabolites"] <- "Detergent"

orderColsBy <- graphData %>%
  group_by(category) %>%
  summarise(median = median(meanEAR[meanEAR != 0])) %>%
  arrange(median)

orderedLevels <- orderColsBy$category

if(any(is.na(orderColsBy$median))){
  orderedLevels <- c(orderColsBy$category[is.na(orderColsBy$median)],
                      orderColsBy$category[!is.na(orderColsBy$median)])
}

orderNames <- names(table(select(endPointInfo, intended_target_family)))

orderedLevels <- c(orderNames[!(orderNames %in% orderedLevels)],orderedLevels)

graphData$category <- factor(as.character(graphData$category), levels=orderedLevels)

countNonZero <- graphData %>%
  group_by(category) %>%
  summarise(nonZero = as.character(length(unique(site[meanEAR>0])))) %>%
  data.frame()

bioPlot <- ggplot(graphData)+
  scale_y_log10("Maximum EAR Per Site",labels=fancyNumbers)+
  geom_boxplot(aes(x=category, y=meanEAR),lwd=0.1,outlier.size=1, fill = "steelblue") +
  coord_flip() +
  theme_bw() +
  xlab("") +
  theme(plot.background = element_rect(fill = "transparent",colour = NA),
        axis.text.y = element_text(size=10, color = "black", vjust = 0.2), 
        axis.text.x = element_text(size=10, color = "black", vjust = 0, margin = margin(-0.5,0,0,0)),
        axis.title = element_text(size=10))

xmin <- 10^(ggplot_build(bioPlot)$layout$panel_ranges[[1]]$x.range[1])
xmax <- 10^(ggplot_build(bioPlot)$layout$panel_ranges[[1]]$x.range[2])
ymax <- ggplot_build(bioPlot)$layout$panel_ranges[[1]]$y.range[1]

bioPlot <- bioPlot + 
  geom_text(data=countNonZero, aes(x=category, y=xmin,label=nonZero),size=3) 

bioPlot <- bioPlot + 
  annotation_custom(
    grob = textGrob(label = "# Non-Zero Sites", 
                    gp = gpar(cex = 0.4)),
    ymin = log10(xmin),      # Vertical position of the textGrob
    ymax = log10(xmin),
    xmin = 16.85,  # Note: The grobs are positioned outside the plot area
    xmax = 16.85)

bioPlot <- ggplot_gtable(ggplot_build(bioPlot))
bioPlot$layout$clip[bioPlot$layout$name == "panel"] <- "off"

ggsave(bioPlot, #bg = "transparent",
       filename = "bioPlot2.png", 
       height = 4, width = 5)

graphData$site[graphData$site == "MilwaukeeMouth"] <- "Milwaukee"

graphData <- graphData %>%
  mutate(site = factor(site, levels = sitesOrdered[sitesOrdered %in% siteLimits$shortName]))

graphData <- left_join(graphData, select(siteLimits, shortName, Lake), by=c("site" = "shortName"))
graphData$Lake <- gsub(" River","", graphData$Lake)
graphData$Lake[graphData$site == "StRegis"] <- "Lake Ontario"
graphData$Lake <- factor(gsub("Lake ","",graphData$Lake), c("Superior", "Michigan","Huron","Erie", "Ontario"))

heat <- ggplot(data = graphData) +
  geom_tile(aes(x = site, y=category, fill=meanEAR)) +
  theme_bw() +
  theme(axis.text.x = element_text( angle = 90,vjust=0.5,hjust = 1)) +
  ylab("") +
  xlab("") +
  labs(fill="Maximum EAR") +
  scale_fill_gradient( guide = "legend",
                       trans = 'log',
                       low = "white", high = "steelblue",
                       breaks=c(0.00001,0.0001,0.001,0.01,0.1,1,5),
                       na.value = 'transparent',labels=fancyNumbers2) +
  facet_grid(. ~ Lake,scales="free", space="free") +
  theme(strip.text.y = element_text(angle=0, hjust=0, size=7), 
        strip.text.x = element_text(size = 7),
        strip.background = element_rect(fill="transparent", colour = NA),
        axis.text = element_text(size=7),
        panel.spacing = unit(0.05, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size=8),
        legend.title = element_text(size =8),
        plot.background = element_rect(fill = "transparent",colour = NA))
heat

ggsave(heat, #bg = "transparent",
       filename = "heat_Bio.png", 
       height = 7, width = 11)
