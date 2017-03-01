library(dplyr)
library(tidyr)

source("getDataReady.R")

intended_target <- select(endPointInfo, intended_target_family, intended_target_family_sub, endPoint = assay_component_endpoint_name, source = assay_source_long_name) %>%
  right_join(select(ep, endPoint), by = "endPoint") %>%
  arrange(intended_target_family, intended_target_family_sub) 

intended_target$intended_target_family_sub["Zebrafish" == intended_target$intended_target_family] <- "Zebrafish"

intended_target <- intended_target %>%
  rename(`Intended Target Family`=intended_target_family,
         `Intended Target Family Sub-Family` = intended_target_family_sub) %>%
  data.frame()

write.csv(intended_target, "intended_target.csv", row.names = FALSE)

chemicalSummary_1 <- chemicalSummary %>%
  left_join(select(endPointInfo, 
                   endPoint=assay_component_endpoint_name,
                   intended_target_family_sub), by="endPoint") %>%
  rename(subFamily = intended_target_family_sub) 

tableData <- chemicalSummary_1 %>%
  group_by(site, Family=choices, subFamily, date, category) %>%
  summarize(sumEAR = sum(EAR)) %>%
  group_by(site, Family, subFamily, category) %>%
  summarize(meanEAR = max(sumEAR)) %>%
  # hits = any(hits > 0)) %>% #is a hit when any EAR is greater than 0.1?
  group_by(Family, subFamily, category) %>%
  summarize(nSites = sum(meanEAR > 10^-3)) %>%
  filter(nSites >= 10) %>%
  data.frame() %>%
  spread(category, nSites) %>%
  arrange(Family) %>%
  select(Family, subFamily, `Bisphenol A`, Metolachlor, Atrazine,
         `Tris(2-butoxyethyl) phosphate`, Caffeine, Cotinine,
         Benzophenone, `Diethyl phthalate`, `Triphenyl phosphate`,
         `Benzo(a)pyrene`,`4-Nonylphenol`)

write.csv(tableData, "intendedTargetCounts.csv", row.names = FALSE, na = "")

#####################################################

newGraph <- chemicalSummary_1 %>%
  group_by(site, choices, subFamily, date) %>%
  summarize(sumEAR = sum(EAR)) %>%
  group_by(site, choices, subFamily) %>%
  summarize(meanEAR = max(sumEAR)) %>%
  data.frame() %>%
  mutate(subFamily = as.character(sapply(tolower(subFamily),simpleCap)))

newGraph$subFamily[newGraph$subFamily == "NANA"] <- "Zebra"

orderSub <- newGraph %>%
  group_by(choices) %>%
  summarise(median = median(meanEAR[meanEAR != 0])) %>%
  data.frame() %>%
  arrange(desc(median))

orderGroups <- newGraph %>%
  group_by(subFamily, choices) %>%
  summarise(median = quantile(meanEAR[meanEAR != 0],0.5)) %>%
  data.frame() %>%
  mutate(choices = factor(choices, levels=orderSub$choices)) %>%
  arrange(choices, desc(median))

newGraph$subFamily <- factor(newGraph$subFamily, levels = rev(orderGroups$subFamily))
newGraph$choices <- factor(newGraph$choices, levels = orderSub$choices)

cbValues <- c("#DCDA4B","#999999","#00FFFF","#CEA226","#CC79A7","#4E26CE",
              "#FFFF00","#78C15A","#79AEAE","#FF0000","#00FF00","#B1611D",
              "#FFA500","#F4426e", "#4286f4","red")

subPlot <- ggplot(newGraph)+
  scale_y_log10("Maximum EAR Per Site",labels=fancyNumbers)+
  geom_boxplot(aes(x=subFamily, y=meanEAR, fill = choices),lwd=0.1,outlier.size=1) +
  coord_flip() +
  theme_bw() +
  xlab("") +
  theme(plot.background = element_rect(fill = "transparent",colour = NA),
        axis.text.y = element_text(size=8, color = "black", vjust = 0.2), 
        axis.text.x = element_text(size=8, color = "black", vjust = 0, margin = margin(-0.5,0,0,0)),
        axis.title = element_text(size=10)) +
  scale_fill_manual(values = cbValues, drop=TRUE)  +
  guides(fill=guide_legend(ncol=6)) +
  theme(legend.position="bottom",
        legend.justification = "left",
        legend.background = element_rect(fill = "transparent", colour = "transparent"),
        legend.title=element_blank(),
        legend.text = element_text(size=8),
        legend.key.height = unit(1,"line")) 

subPlot

ggsave(subPlot, bg = "transparent",
       filename = "subFamilies.png", 
       height = 10, width = 7.75)

#####################################################

tableData_EPs <- chemicalSummary %>%
  group_by(site, choices, endPoint, date, category) %>%
  summarize(sumEAR = sum(EAR)) %>%
  group_by(site, choices, endPoint, category) %>%
  summarize(meanEAR = max(sumEAR)) %>%
  # hits = any(hits > 0)) %>% #is a hit when any EAR is greater than 0.1?
  group_by(choices, endPoint, category) %>%
  summarize(nSites = sum(meanEAR > 10^-3)) %>%
  # filter(nSites >= 10) %>%
  data.frame() %>%
  spread(category, nSites) %>%
  arrange(choices) %>%
  select(choices, endPoint, `Bisphenol A`, Metolachlor, Atrazine,
         `Tris(2-butoxyethyl) phosphate`, Caffeine, Cotinine,
         Benzophenone, `Diethyl phthalate`, `Triphenyl phosphate`,
         `4-Nonylphenol`)

write.csv(tableData, "intendedTarget.csv", row.names = FALSE, na = "")


tableData <- chemicalSummary_1 %>%
  group_by(site, endPoint, Family=choices, subFamily, date, Chemical=category) %>%
  summarize(sumEAR = sum(EAR)) %>%
  group_by(site, endPoint, Family, subFamily, Chemical) %>%
  summarize(meanEAR = max(sumEAR)) %>%
  # hits = any(hits > 0)) %>% #is a hit when any EAR is greater than 0.1?
  group_by(endPoint, Family, subFamily, Chemical) %>%
  summarize(nSites = sum(meanEAR > 10^-3)) %>%
  # filter(nSites >= 10) %>%
  data.frame() %>%
  filter(nSites > 0) %>%
  spread(Chemical, nSites) %>%
  arrange(Family, subFamily) %>%
  select(endPoint, Family, subFamily, everything())

sortedCols <- sort(colSums(tableData[,4:length(tableData)], na.rm = TRUE), 
     decreasing = TRUE)
tableData <- tableData[,c("endPoint","Family","subFamily", names(sortedCols))]

tableData2 <- select(tableData, -endPoint, -Family, -subFamily) 

tableData$nChems <- apply(tableData2, MARGIN = 1, function(x) sum(x>0, na.rm = TRUE))

tableData3 <- tableData[,c("Family", "subFamily", "endPoint", "nChems", rev(orderChem$category)[rev(orderChem$category) %in% names(tableData)])]

write.csv(tableData3, row.names = FALSE, file="wholeThing.csv", na = "")

