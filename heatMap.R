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

graphData$class <- factor(graphData$class, levels=c(levels(orderClass$class),""))
graphData$category <- factor(graphData$category, levels=c("Total","Class Summation",orderedLevels))

siteLimits <- siteLimits %>%
  mutate(shortName = factor(shortName, levels=sitesOrdered[sitesOrdered %in% siteLimits$shortName]))

graphData <- left_join(graphData, select(siteLimits, shortName, Lake), by=c("site" = "shortName"))
graphData$Lake <- gsub(" River","", graphData$Lake)
graphData$Lake <- factor(gsub("Lake ","",graphData$Lake), c("Superior", "Michigan","Huron","Erie", "Ontario"))

sitesToUse <- graphData %>%
  group_by(site) %>%
  summarize(toUse = any(meanEAR > 10^-3))

classSum <- chemicalSummary %>%
  group_by(site,date,category,class) %>% #sample
  summarise(sumEAR=sum(EAR)) %>%
  data.frame() %>%
  group_by(site, class, category) %>% #max at site
  summarise(maxEAR=max(sumEAR)) %>%
  group_by(site, class) %>%
  summarise(meanEAR=sum(maxEAR)) %>%
  data.frame() %>%
  mutate(category = factor("Class Summation", levels = c("Total","Class Summation",orderedLevels)),
         site = factor(site, levels = sitesOrdered[sitesOrdered %in% siteLimits$shortName]),
         class = factor(class, levels = c("",levels(orderClass$class))))

overallSum <- chemicalSummary %>%
  group_by(site,date) %>% 
  summarise(sumEAR=sum(EAR)) %>%
  data.frame() %>%
  group_by(site) %>% #max at site
  summarise(meanEAR=max(sumEAR)) %>%
  data.frame() %>%
  mutate(category = factor("Total", levels = c("Total","Class Summation",orderedLevels)),
         site = factor(site, levels = sitesOrdered[sitesOrdered %in% siteLimits$shortName]),
         class = factor("", levels = c(levels(orderClass$class),"")))

classSum <- left_join(classSum, select(siteLimits, shortName, Lake), by=c("site" = "shortName"))
classSum$Lake <- gsub(" River","", classSum$Lake)
classSum$Lake <- factor(gsub("Lake ","",classSum$Lake), c("Superior", "Michigan","Huron","Erie", "Ontario"))

overallSum <- left_join(overallSum, select(siteLimits, shortName, Lake), by=c("site" = "shortName"))
overallSum$Lake <- gsub(" River","", overallSum$Lake)
overallSum$Lake <- factor(gsub("Lake ","",overallSum$Lake), c("Superior", "Michigan","Huron","Erie", "Ontario"))

graphData <- bind_rows(graphData, classSum, overallSum)

graphData_filtered <- graphData %>%
  filter(site %in% sitesToUse$site[sitesToUse$toUse]) 

classToUse <- filter(graphData_filtered, !(category %in% c("Total","Class Summation"))) %>%
  group_by(class, site, category) %>%
  summarize(toUse = any(meanEAR > 10^-3)) %>%
  group_by(class) %>%
  summarize(with10 = sum(toUse) >= 20,
            countCats = length(unique(category))) %>%
  filter(countCats > 2)

graphData_filtered <- graphData_filtered %>%
  filter(class %in% classToUse$class[classToUse$with10])

overallSum_filtered <- overallSum %>%
  filter(site %in% sitesToUse$site[sitesToUse$toUse]) 

graphData_filtered <- bind_rows(graphData_filtered, overallSum_filtered)

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
  facet_grid(class ~ Lake,scales="free", space="free") +
  theme(strip.text.y = element_text(angle=0, hjust=0, size=7), 
        strip.text.x = element_text(size = 8),
        strip.background = element_rect(fill="transparent", colour = NA),
        axis.text = element_text(size=9),
        axis.text.y = element_text(face=ifelse(levels(graphData_filtered$category) %in% c("Total"),"bold","italic")),
        panel.spacing = unit(0.05, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size=8),
        legend.title = element_text(size =8),
        plot.background = element_rect(fill = "transparent",colour = NA))
heat

ggsave(heat, #bg = "transparent",
       filename = "heat_2.png", 
       height = 5, width = 9)

ggsave(heat, #bg = "transparent",
       filename = "heat_Unfiltered.png", 
       height = 9, width = 11)

