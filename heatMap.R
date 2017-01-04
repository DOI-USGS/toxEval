library(toxEval)
library(dplyr)
library(data.table)
library(tidyr)
library(ggplot2)
library(grid)
library(gridExtra)

source("getDataReady.R")

graphData <- graphData %>%
  mutate(site = factor(site, levels = sitesOrdered[sitesOrdered %in% siteLimits$shortName])) %>%
  filter(!is.na(category)) 

graphData$class <- factor(graphData$class, levels=levels(orderClass$class))
graphData$category <- factor(graphData$category, levels=orderedLevels)

siteLimits <- siteLimits %>%
  mutate(shortName = factor(shortName, levels=sitesOrdered[sitesOrdered %in% siteLimits$shortName]))

graphData <- left_join(graphData, select(siteLimits, shortName, Lake), by=c("site" = "shortName"))
graphData$Lake <- gsub(" River","", graphData$Lake)
graphData$Lake <- factor(gsub("Lake ","",graphData$Lake), c("Superior", "Michigan","Huron","Erie", "Ontario"))

sitesToUse <- graphData %>%
  group_by(site) %>%
  summarize(toUse = any(meanEAR > 10^-3))

graphData_filtered <- graphData %>%
  filter(site %in% sitesToUse$site[sitesToUse$toUse]) 

classToUse <- graphData_filtered %>%
  group_by(class, site, category) %>%
  summarize(toUse = any(meanEAR > 10^-3)) %>%
  group_by(class) %>%
  summarize(with10 = sum(toUse) >= 20,
            countCats = length(unique(category))) %>%
  filter(countCats > 2)

graphData_filtered <- graphData_filtered %>%
  filter(class %in% classToUse$class[classToUse$with10])

# graphData_filtered <- graphData %>%
#   filter(!(Lake %in% c("Superior","Huron","Ontario","St. Lawrence"))) %>%
#   filter(!is.na(category)) %>%
#   filter(class %in% c("Plasticizer", "Herbicide","PAH", "Antioxidant"))

# rmCat <- group_by(graphData_filtered, category) %>%
#   summarize(allZero = all(meanEAR == 0))

heat <- ggplot(data = graphData_filtered) +
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
  facet_grid(class ~ Lake,scales="free", space="free") +
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

ggsave(heat, bg = "transparent",
       filename = "heat_2.png", 
       height = 5, width = 8.5)

ggsave(heat, bg = "transparent",
       filename = "heat_Unfiltered.png", 
       height = 9, width = 11)

