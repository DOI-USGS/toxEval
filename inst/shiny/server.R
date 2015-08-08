library(dplyr)
library(magrittr)
library(ggplot2)
library(DT)
library(leaflet)
library(toxEval)
library(data.table)
library(gridExtra)
library(grid)
library(tidyr)
library(RColorBrewer)

endPointInfo <- endPointInfo

packagePath <- system.file("extdata", package="toxEval")
filePath <- file.path(packagePath, "stationINFO.RData")
load(file=filePath)
siteKey <- setNames(stationINFO$shortName, stationINFO$fullSiteID)
stationINFO <- stationINFO %>%
  mutate(lakeCat = factor(Lake, 
                          levels=c("Lake Superior","Lake Michigan",
                                   "Lake Huron",
                                   "Detroit River and Lake St. Clair","Lake Erie",
                                   "Lake Ontario", "St. Lawrence River"))) %>%
  arrange(lakeCat, dec.long.va) %>%
  mutate(lakeColor = c("tomato3","black","springgreen3","brown","brown","blue","blue")[as.numeric(lakeCat)] )

lakeKey <- setNames(as.character(stationINFO$lakeCat), stationINFO$fullSiteID)
lakeKey[lakeKey == "Detroit River and Lake St. Clair"] <- "Lake Erie"
lakeKey[lakeKey == "St. Lawrence River"] <- "Lake Ontario"


endPointInfo <- endPointInfo

pathToApp <- system.file("extdata", package="toxEval")

summaryFile <- readRDS(file.path(pathToApp,"summary.rds"))
endPoint <- readRDS(file.path(pathToApp,"endPoint.rds"))
chemicalSummaryWS <- readRDS(file.path(pathToApp,"chemicalSummary.rds"))
chemicalSummaryPS <- readRDS(file.path(pathToApp,"chemicalSummaryPassive.rds"))
chemicalSummaryPS$date <- rep(as.POSIXct(as.Date("1970-01-01")),nrow(chemicalSummaryPS))
df2016 <- readRDS(file.path(pathToApp,"df2016.rds"))

choicesPerGroup <- apply(endPointInfo[,-3], 2, function(x) length(unique(x)))
groupChoices <- paste0(names(choicesPerGroup)," (",choicesPerGroup,")")

interl3 <- function (a,b,cv) {
  n <- min(length(a),length(b),length(cv))
  p1 <- as.vector(rbind(a[1:n],b[1:n],cv[1:n]))
  p2 <- c(a[-(1:n)],b[-(1:n)],cv[-(1:n)])
  c(p1,p2)
}

interl <- function (a,b) {
  n <- min(length(a),length(b))
  p1 <- as.vector(rbind(a[1:n],b[1:n]))
  p2 <- c(a[-(1:n)],b[-(1:n)])
  c(p1,p2)
}

shinyServer(function(input, output) {
  
#############################################################   
    chemicalSummary <- reactive({
      
      if(is.null(input$data)){
        chemicalSummary <- rbind(chemicalSummaryWS, chemicalSummaryPS)
      } else if (input$data == "Passive Samples"){
        chemicalSummary <- chemicalSummaryPS
      } else if (input$data == "Water Sample"){
        chemicalSummary <- chemicalSummaryWS
      } else {
        chemicalSummary <- rbind(chemicalSummaryWS, chemicalSummaryPS)
      }
      
      if(is.null(input$sites) | input$sites == "All"){
        siteToFind <- summaryFile$site
      } else if (input$sites == "Potential 2016"){
        siteToFind <- df2016$Site
      } else {
        siteToFind <- input$sites
      }
      
      chemicalSummary 
      
    })
  
    chemicalSummaryFiltered <- reactive({
      
      chemicalSummary <- chemicalSummary()
      
      if(is.null(input$groupCol)){
        groupCol <- names(endPointInfo)[20]
      } else {
        groupCol <- input$groupCol
      }

      endPointInfoSub <- select_(endPointInfo, "assay_component_endpoint_name", groupCol) %>%
        distinct()
      
      chemicalSummaryFiltered <- chemicalSummary %>%
        rename(assay_component_endpoint_name=endPoint) %>%
        filter(assay_component_endpoint_name %in% endPointInfoSub$assay_component_endpoint_name ) %>%
        data.table()%>%
        left_join(data.table(endPointInfoSub), by = "assay_component_endpoint_name") %>%
        data.frame()%>%
        rename(endPoint=assay_component_endpoint_name)

      if(nrow(chemicalSummaryFiltered) == 0){
        endPointInfoSub <- select_(endPointInfo, 
                                   "assay_component_endpoint_name", names(endPointInfo)[20]) %>%
          distinct()
        chemicalSummaryFiltered <- chemicalSummary() %>%
          rename(assay_component_endpoint_name=endPoint) %>%
          filter(assay_component_endpoint_name %in% endPointInfoSub$assay_component_endpoint_name ) %>%
          data.table()%>%
          left_join(data.table(endPointInfoSub), by = "assay_component_endpoint_name") %>%
          data.frame()%>%
          rename(endPoint=assay_component_endpoint_name)
        
      }      
      
      chemicalSummaryFiltered
      
    })

    endpointSummary <- reactive({
      
      chemicalSummaryFiltered <- chemicalSummaryFiltered()
      
      if(is.null(input$groupCol)){
        groupCol <- names(endPointInfo)[20]
      } else {
        groupCol <- input$groupCol
      }
      
      if(is.null(input$group)){
        group <- unique(endPointInfo[,20])[3]
      } else {
        group <- input$group
      }
      
      endpointSummary <- chemicalSummaryFiltered %>%
        filter_(paste0(groupCol," == '", group, "'")) %>%
        group_by(site, date) %>%
        summarise(sumEAR = sum(EAR),
                  nHits = sum(hits)) 
    })
    
    statsOfColumn <- reactive({
      
      chemicalSummary <- chemicalSummary()
      
      if(is.null(input$sites) | input$sites == "All"){
        siteToFind <- summaryFile$site
      } else if (input$sites == "Potential 2016"){
        siteToFind <- df2016$Site
      } else {
        siteToFind <- input$sites
      }
      
      if(is.null(input$groupCol)){
        groupCol <- names(endPointInfo)[20]
      } else {
        groupCol <- input$groupCol
      }
      
      if(is.null(input$radioMaxGroup)){
        radio <- "2"
      } else {
        radio <- input$radioMaxGroup
      }
      
      statsOfColumn <- chemicalSummary %>%
        rename(assay_component_endpoint_name=endPoint) %>%
        filter(assay_component_endpoint_name %in% endPointInfo$assay_component_endpoint_name ) %>%
        data.table()%>%
        left_join(data.table(endPointInfo[,c("assay_component_endpoint_name", groupCol)]), by = "assay_component_endpoint_name") %>%
        data.frame() %>%
        mutate(site = siteKey[site]) %>%
        select_("hits","EAR","chnm","class","date","choices"=groupCol,"site") %>%
        filter(site %in% siteToFind)
      
      if(length(siteToFind) == 1){
 
        statsOfColumn <- statsOfColumn %>%
          # mutate(site=siteKey[site]) %>%
          filter(site %in% siteToFind) 
        
        if(radio == "2"){
          statsOfColumn <- mutate(statsOfColumn, site = chnm) 
        } else if (radio == "3"){
          statsOfColumn <- mutate(statsOfColumn, site = class) 
        } else {
          statsOfColumn <- mutate(statsOfColumn, site = choices) 
        }
      } 

      statsOfColumn <- statsOfColumn %>%
        group_by(site, date, choices) %>%
        summarise(sumEAR = sum(EAR),
                  nHits = sum(hits)) %>%
        group_by(site, choices) %>%
        summarise(maxEAR = max(sumEAR),
                  sumHits = sum(nHits),
                  freq = sum(nHits > 0)/n()) %>%
        data.frame()
      
      if(!(length(siteToFind) == 1 & radio == "1")){
        statsOfColumn <- statsOfColumn %>%
          gather(calc, value, -site, -choices) %>%
          unite(choice_calc, choices, calc, sep=" ") %>%
          spread(choice_calc, value)        
      }
      
      if(length(siteToFind) > 1){
        if (input$data == "Water Sample"){
          summaryFile <- filter(summaryFile, site %in% siteToFind)
          
          statsOfColumn <- left_join(summaryFile[,c("site","nSamples")], statsOfColumn, by="site")
        }
      }
      
      statsOfColumn
      
    })
    
    chemGroup <- reactive({
      
      chemicalSummary <- chemicalSummary()
      
      if(is.null(input$groupCol)){
        groupCol <- names(endPointInfo)[20]
      } else {
        groupCol <- input$groupCol
      }

      chemGroup <- chemicalSummary %>%
        rename(assay_component_endpoint_name=endPoint) %>%
        filter(assay_component_endpoint_name %in% endPointInfo$assay_component_endpoint_name ) %>%
        data.table()%>%
        left_join(data.table(endPointInfo[,c("assay_component_endpoint_name", groupCol)]), by = "assay_component_endpoint_name") %>%
        data.frame()%>%
        mutate(site = siteKey[site]) %>%
        rename(endPoint=assay_component_endpoint_name) %>%
        select_("hits","EAR","chnm","class","site","date","choices"=groupCol)
      
    })
    
    statsOfGroup <- reactive({
      
      if(is.null(input$groupCol)){
        groupCol <- names(endPointInfo)[20]
      } else {
        groupCol <- input$groupCol
      }
      
      if(is.null(input$sites) | input$sites == "All"){
        siteToFind <- summaryFile$site
      } else if (input$sites == "Potential 2016"){
          siteToFind <- df2016$Site
      } else {
        siteToFind <- input$sites
      }
      
      if(is.null(input$radioMaxGroup)){
        radioMaxGroup <- "1"
      } else {
        radioMaxGroup <- input$radioMaxGroup
      }
      
      chemGroup <- chemGroup()   
      
      statsOfGroup <-  chemGroup %>%
        filter(site %in% siteToFind)
      
      # Move this:      
      if(radioMaxGroup == "2"){
        statsOfGroup <- statsOfGroup %>%
          mutate(category=chnm)
      } else if(radioMaxGroup == "3"){
        statsOfGroup <- statsOfGroup %>%
          mutate(category=class)        
      } else {
        statsOfGroup <- statsOfGroup %>%
          mutate(category = choices)         
      }

      statsOfGroup <- statsOfGroup %>%
        group_by(site, date,category) %>%
        summarise(hits = length(category[EAR > 0.1]),
                  sumEAR = sum(EAR)) %>%
        group_by(site,category) %>%
        summarise(count = sum(hits > 0),
                  max = max(sumEAR),
                  mean = mean(sumEAR))
      
      statsOfGroup
        

    })
    
    statsOfGroupOrdered <- reactive({
      
      statsOfGroup <- statsOfGroup() 
      
      if(is.null(input$sites) | input$sites == "All"){
        siteToFind <- summaryFile$site
      } else if (input$sites == "Potential 2016"){
        siteToFind <- df2016$Site
      } else {
        siteToFind <- input$sites
      }
      
      statsOfGroup <- statsOfGroup %>%
        filter(site %in% siteToFind)
      
      if(length(siteToFind) > 1){
        
        statsOfGroupOrdered <- statsOfGroup %>%
          gather(calc, value, -site, -category) %>%
          unite(choice_calc, category, calc, sep=" ") %>%
          spread(choice_calc, value)
        
        maxCat <- grep("max",names(statsOfGroupOrdered))
        
        maxCatordered <- order(apply(statsOfGroupOrdered[,maxCat], 2, max),decreasing = TRUE)
        
        if(length(maxCat) > 9){
          statsOfGroupOrdered <- statsOfGroupOrdered[,c(1,interl3((maxCat[maxCatordered[1:9]]-1),maxCat[maxCatordered[1:9]],(maxCat[maxCatordered[1:9]]+1)))]
          maxCat <- maxCat[1:9]
        } else {
          if(length(maxCat) > 1){
            statsOfGroupOrdered <- statsOfGroupOrdered[,c(1,interl3((maxCat[maxCatordered]-1),maxCat[maxCatordered],(maxCat[maxCatordered]+1)))]
          }
        }
        
      } else {
        statsOfGroupOrdered <- statsOfGroup
      }
      
      statsOfGroupOrdered
        
    }) 
    
    chemGroupBP_group <- reactive({
      
      if(is.null(input$radioMaxGroup)){
        radioMaxGroup <- "1"
      } else {
        radioMaxGroup <- input$radioMaxGroup
      }
      
      chemGroup <- chemGroup()
      
      chemGroupBP_group <- chemGroup %>%
        # data.frame() %>%
        mutate(date = factor(as.numeric(date) - min(as.numeric(date))))
      
      if(radioMaxGroup == "1"){
        chemGroupBP_group <- mutate(chemGroupBP_group, cat=choices)
      } else if (radioMaxGroup == "2"){
        chemGroupBP_group <- mutate(chemGroupBP_group, cat=chnm)
      } else {
        chemGroupBP_group <- mutate(chemGroupBP_group, cat=class)
      }
      
      chemGroupBP_group     
      
    })
    
    chemGroupBP <- reactive({
      
      if(is.null(input$radio)){
        radioMaxGroup <- "1"
      } else {
        radioMaxGroup <- input$radio
      }
      
      chemGroup <- chemGroup()
      
      chemGroupBP <- chemGroup %>%
        # data.frame() %>%
        mutate(date = factor(as.numeric(date) - min(as.numeric(date))))
      
      if(radioMaxGroup == "2"){
        chemGroupBP <- mutate(chemGroupBP, cat=class)
      } else {
        chemGroupBP <- mutate(chemGroupBP, cat=chnm)
      } 
      
      chemGroupBP     
      
    })
    
############################################################# 
        
    output$table <- DT::renderDataTable({

      if(is.null(input$radio)){
        radio <- 1
      } else {
        radio <- input$radio
      }
      
      if(is.null(input$groupCol)){
        groupCol <- names(endPointInfo)[20]
      } else {
        groupCol <- input$groupCol
      }
      
      if(is.null(input$group)){
        group <- unique(endPointInfo[,20])[3]
      } else {
        group <- input$group
      }
      
      if(is.null(input$sites) | input$sites == "All"){
        siteToFind <- summaryFile$site
      } else if (input$sites == "Potential 2016"){
        siteToFind <- df2016$Site
      } else {
        siteToFind <- input$sites
      }
      
      chemicalSummaryFiltered <- chemicalSummaryFiltered()
      
      if(radio == "2"){
        chemicalSummaryFiltered <- chemicalSummaryFiltered %>%
          mutate(grouping=class)
      } else {
        chemicalSummaryFiltered <- chemicalSummaryFiltered %>%
          mutate(grouping=chnm)        
      }
      
      statTable <- chemicalSummaryFiltered %>%
        filter_(paste0(groupCol," == '", group, "'")) 
        mutate(site = siteKey[site]) %>%
        filter(site %in% siteToFind) %>%
        group_by(grouping, date) %>%
        summarise(sumEAR = sum(EAR),
                  nHits = sum(hits)) %>%
        group_by(grouping) %>%
        summarise(meanEAR = mean(sumEAR),
                  maxEAR = max(sumEAR),
                  freq = sum(nHits > 0)/n()) %>%
        data.frame()

        statsOfSumDF <- DT::datatable(statTable, 
                                      rownames = FALSE,
                                      colnames = c('Maximum EAR' = 3, 'Frequency of Hits' = 4,
                                                   'Mean EAR' = 2),
                                      filter = 'top',
                                      options = list(pageLength = nrow(statTable), 
                                                     order=list(list(2,'desc')))) %>%
          formatRound(c("Maximum EAR","Frequency of Hits","Mean EAR"), 2) %>%
          formatStyle("Maximum EAR", 
                      background = styleColorBar(range(statTable[,3],na.rm = TRUE), 'goldenrod'),
                      backgroundSize = '100% 90%',
                      backgroundRepeat = 'no-repeat',
                      backgroundPosition = 'center'
          ) %>%
          formatStyle("Frequency of Hits", 
                      background = styleColorBar(range(statTable[,4], na.rm=TRUE), 'wheat'),
                      backgroundSize = '100% 90%',
                      backgroundRepeat = 'no-repeat',
                      backgroundPosition = 'center'
          ) %>%
          formatStyle("Mean EAR", 
                      background = styleColorBar(range(statTable[,2], na.rm=TRUE), 'seashell'),
                      backgroundSize = '100% 90%',
                      backgroundRepeat = 'no-repeat',
                      backgroundPosition = 'center'
          )

    })

    output$tableGroupSumm <- DT::renderDataTable({        

      statsOfGroupOrdered <- statsOfGroupOrdered()
  
      if(is.null(input$sites) | input$sites == "All"){
        siteToFind <- summaryFile$site
      } else if (input$sites == "Potential 2016"){
        siteToFind <- df2016$Site
      } else {
        siteToFind <- input$sites
      }

      statsOfGroupOrdered <- statsOfGroupOrdered %>%
        filter(site %in% siteToFind)
      
      if(length(siteToFind) > 1){
        
        if (input$data == "Water Sample"){
          
          summaryFile <- filter(summaryFile, site %in% siteToFind)
          
          statsOfGroupOrdered <- left_join(summaryFile[,c("site","nSamples")], statsOfGroupOrdered, by="site")
        }
        
        meanChem <- grep("mean",names(statsOfGroupOrdered))
        maxChem <- grep("max",names(statsOfGroupOrdered))
        nChem <- grep(" count",names(statsOfGroupOrdered))
        
        colors <- brewer.pal(length(maxChem),"Blues") #"RdYlBu"

        tableGroup <- DT::datatable(statsOfGroupOrdered, 
                                      rownames = FALSE,
                                      filter = 'top',
                                      options = list(pageLength = nrow(statsOfGroupOrdered), 
                                                     order=list(list(2,'desc'))))
        tableGroup <- formatRound(tableGroup, names(statsOfGroupOrdered)[c(meanChem, maxChem)], 2)
          
          
          for(i in 1:length(maxChem)){
            tableGroup <- formatStyle(tableGroup, 
                                      names(statsOfGroupOrdered)[maxChem[i]], 
                                      backgroundColor = colors[i])
            tableGroup <- formatStyle(tableGroup, 
                                      names(statsOfGroupOrdered)[meanChem[i]], 
                                      backgroundColor = colors[i])
            tableGroup <- formatStyle(tableGroup, 
                                      names(statsOfGroupOrdered)[nChem[i]], 
                                      backgroundColor = colors[i])
            
            tableGroup <- formatStyle(tableGroup, names(statsOfGroupOrdered)[maxChem[i]], 
                                      background = styleColorBar(range(statsOfGroupOrdered[,maxChem[i]],na.rm=TRUE), 'goldenrod'),
                                      backgroundSize = '100% 90%',
                                      backgroundRepeat = 'no-repeat',
                                      backgroundPosition = 'center' ) 
            
            tableGroup <- formatStyle(tableGroup, names(statsOfGroupOrdered)[meanChem[i]], 
                                      background = styleColorBar(range(statsOfGroupOrdered[,meanChem[i]],na.rm=TRUE), 'wheat'),
                                      backgroundSize = '100% 90%',
                                      backgroundRepeat = 'no-repeat',
                                      backgroundPosition = 'center') 
            
            tableGroup <- formatStyle(tableGroup, names(statsOfGroupOrdered)[nChem[i]], 
                                      background = styleColorBar(range(statsOfGroupOrdered[,nChem[i]],na.rm=TRUE), 'seashell'),
                                      backgroundSize = '100% 90%',
                                      backgroundRepeat = 'no-repeat',
                                      backgroundPosition = 'center') 
          }
          
        } else {
          
          statsOfGroupOrdered <- statsOfGroupOrdered[, -1]
          
          tableGroup <- DT::datatable(statsOfGroupOrdered, 
                                      rownames = FALSE,
                                      filter = 'top',
                                      options = list(pageLength = nrow(statsOfGroupOrdered), 
                                                     order=list(list(1,'desc'))))
          tableGroup <- formatRound(tableGroup, c("max","mean"), 2)
          tableGroup <- formatStyle(tableGroup, "max", 
                                    background = styleColorBar(range(statsOfGroupOrdered[,"max"],na.rm=TRUE), 'goldenrod'),
                                    backgroundSize = '100% 90%',
                                    backgroundRepeat = 'no-repeat',
                                    backgroundPosition = 'center' ) 
          
          tableGroup <- formatStyle(tableGroup, "mean", 
                                    background = styleColorBar(range(statsOfGroupOrdered[,"mean"],na.rm=TRUE), 'wheat'),
                                    backgroundSize = '100% 90%',
                                    backgroundRepeat = 'no-repeat',
                                    backgroundPosition = 'center') 
        }
        
        tableGroup
      
    })
    
    output$tableSumm <- DT::renderDataTable({
      
      statCol <- statsOfColumn()
      
      freqCol <- grep("freq",names(statCol))
      maxEARS <- grep("maxEAR",names(statCol))
      
      ignoreIndex <- which(names(statCol) %in% c("site","nSamples"))
      
      statCol <- statCol[,c(ignoreIndex,c(maxEARS,freqCol)[order(c(maxEARS,freqCol))])]
      
      maxEARS <- grep("maxEAR",names(statCol))
      
      MaxEARSordered <- order(apply(statCol[,maxEARS, drop = FALSE], 2, max),decreasing = TRUE)
      
      if(length(maxEARS) > 9){
        statCol <- statCol[,c(ignoreIndex,interl(maxEARS[MaxEARSordered[1:9]],(maxEARS[MaxEARSordered[1:9]]-1)))]
        maxEARS <- maxEARS[1:9]
      } else {
        if(length(maxEARS) != 1){
          statCol <- statCol[,c(ignoreIndex,interl(maxEARS[MaxEARSordered],(maxEARS[MaxEARSordered]-1)))]
        }
      }
      
      freqCol <- grep("freq",names(statCol))
      maxEARS <- grep("maxEAR",names(statCol))
      
      colors <- brewer.pal(length(maxEARS),"Blues") #"RdYlBu"
      tableSumm <- DT::datatable(statCol, 
                                 rownames = FALSE,
                                 options = list(pageLength = nrow(statCol),
                                   order=list(list(2,'desc'))))

      
      tableSumm <- formatRound(tableSumm, names(statCol)[-ignoreIndex], 2) 
      
      for(i in 1:length(maxEARS)){
        tableSumm <- formatStyle(tableSumm, 
                                 names(statCol)[maxEARS[i]], 
                                 backgroundColor = colors[i])
        tableSumm <- formatStyle(tableSumm, 
                                 names(statCol)[freqCol[i]], 
                                 backgroundColor = colors[i])
        
        tableSumm <- formatStyle(tableSumm, names(statCol)[maxEARS[i]], 
                                 background = styleColorBar(range(statCol[,names(statCol)[maxEARS[i]]],na.rm = TRUE), 'goldenrod'),
                                 backgroundSize = '100% 90%',
                                 backgroundRepeat = 'no-repeat',
                                 backgroundPosition = 'center' ) 
        tableSumm <- formatStyle(tableSumm, names(statCol)[freqCol[i]], 
                                 background = styleColorBar(range(statCol[,names(statCol)[freqCol[i]]],na.rm = TRUE), 'wheat'),
                                 backgroundSize = '100% 90%',
                                 backgroundRepeat = 'no-repeat',
                                 backgroundPosition = 'center') 
        
      }
      
      tableSumm
      
    })
    
#############################################################
    
    output$stackBarGroup <- renderPlot({ 
      print(topPlotsGroup())
    })

    output$stackBar <- renderPlot({ 
      print(topPlots())
    })

    output$graph <- renderPlot({ 
      print(bottomPlots())
    })
    
    output$graphGroup <- renderPlot({ 
      print(bottomPlotsGroups())
    })
    
    bottomPlots <- reactive({
      
      if(is.null(input$sites) | input$sites == "All"){
        siteToFind <- summaryFile$site
      } else if (input$sites == "Potential 2016"){
        siteToFind <- df2016$Site
      } else {
        siteToFind <- input$sites
      }
      
      if(is.null(input$radioMaxGroup)){
        radio <- 2
      } else {
        radio <- input$radio
      }
      
      if(length(siteToFind) > 1){
        
        endpointSummary <- endpointSummary()
        
        siteLimits <- stationINFO %>%
          filter(shortName %in% siteToFind) %>%
          filter(fullSiteID %in% unique(endpointSummary$site))
        
        if(nrow(endpointSummary) == 0){
          endpointSummary$site <- stationINFO$fullSiteID
        }
        
        endPointSummBP <- endpointSummary %>%
          filter(site %in% siteToFind) %>%
          data.frame()%>%
          mutate(lake = as.character(lakeKey[site])) %>%
          mutate(lake = factor(lake, levels=c("Lake Superior","Lake Michigan",
                                              "Lake Huron","Lake Erie","Lake Ontario"))) %>%
          mutate(site = siteKey[site]) %>%
          mutate(site = factor(site, levels=siteLimits$Station.shortname)) %>%
          mutate(sumEARnoZero = sumEAR) 
        
        
        ndLevel <- 0.1*min(endPointSummBP$sumEARnoZero[endPointSummBP$sumEARnoZero != 0])
        
        endPointSummBP$sumEARnoZero[endPointSummBP$sumEARnoZero == 0] <- ndLevel
        
        
        sToxWS <- ggplot(endPointSummBP)
        
        if(!is.null(input$data) && input$data != "Passive Samples"){
          sToxWS <- sToxWS + geom_boxplot(aes(x=site, y=sumEARnoZero, fill = lake)) +
            scale_fill_manual(values=c("tomato3",
                                       "black",
                                       "springgreen3",
                                       "brown",
                                       "blue"))
        } else {
          sToxWS <- sToxWS + geom_point(aes(x=site, y=sumEARnoZero, 
                                            colour = lake))+
            scale_colour_manual(values=c("tomato3",
                                         "black",
                                         "springgreen3",
                                         "brown",
                                         "blue"))
        }
        
        sToxWS <- sToxWS + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.25, 
                                                            colour=siteLimits$lakeColor))+
          scale_x_discrete(limits=siteLimits$Station.shortname) +
          scale_y_log10("Summation of EAR per sample") +
          coord_cartesian(ylim = c(0.1, 1.1*max(endPointSummBP$sumEARnoZero))) +
          xlab("") 
        
      } else {
        
        chemicalSummaryFiltered <- chemicalSummaryFiltered()
        
        g <- ggplot_build(topPlots())
        fillColors <- g$data[[1]]["fill"][[1]]
        groupings <- g$plot$data$grouping
        
        dfNames <- data.frame(grouping = names(table(groupings)[table(groupings) != 0]),
                              freq = as.numeric(table(groupings)[table(groupings) != 0]),stringsAsFactors = FALSE)
        dfCol <- data.frame(colors = names(table(fillColors)), 
                            freq = as.numeric(table(fillColors)),stringsAsFactors = FALSE)
        
        colorKey <- left_join(dfNames, dfCol, by="freq")
        
        if(radio == "3"){
          chemicalSummaryFiltered <- chemicalSummaryFiltered %>%
            mutate(grouping=class)
        } else {
          chemicalSummaryFiltered <- chemicalSummaryFiltered %>%
            mutate(grouping=chnm)        
        }
        
        chemSiteSumm <- chemicalSummaryFiltered %>%
          filter_(paste0(groupCol," == '", group, "'")) %>%
          mutate(site = siteKey[site]) %>%
          filter(site %in% siteToFind) %>%
          group_by(grouping, date) %>%
          summarise(sumEAR = sum(EAR),
                    nHits = sum(hits))%>%
          mutate(sumEARnoZero = sumEAR) 
        
        ndLevel <- 0.1*min(chemSiteSumm$sumEARnoZero[chemSiteSumm$sumEARnoZero != 0])
        
        if(is.finite(ndLevel)){
          chemSiteSumm$sumEARnoZero[chemSiteSumm$sumEARnoZero == 0] <- ndLevel
          
          chemicals <- unique(chemSiteSumm$grouping)[!(unique(chemSiteSumm$grouping) %in% colorKey$grouping)]
          
          noFill <- data.frame(grouping=chemicals,
                               colors=rep("#FFFFFF",length(chemicals)),
                               stringsAsFactors = FALSE)
          fillTotal <- rbind(colorKey[,c('colors','grouping')], noFill)
          
          sToxWS <- ggplot(chemSiteSumm) 
          
          if(!is.null(input$data) && input$data != "Passive Samples"){
            sToxWS <- sToxWS + geom_boxplot(aes(x=grouping, y=sumEAR, 
                                                fill=grouping)) +
              scale_fill_manual(values=setNames(fillTotal$colors,fillTotal$grouping))
          } else {
            sToxWS <- sToxWS + geom_point(aes(x=grouping, y=sumEAR)) 
          }
          sToxWS <- sToxWS + theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.25), 
                                   legend.position = "none")+
            xlab("") +
            scale_y_log10("Summation of EAR per sample") +
            geom_hline(yintercept=0.1)
        }
      }
      
      sToxWS
    })

    topPlotsGroup <- reactive({
  
      if(is.null(input$sites) | input$sites == "All"){
        siteToFind <- summaryFile$site
      } else if (input$sites == "Potential 2016"){
        siteToFind <- df2016$Site
      } else {
        siteToFind <- input$sites
      }
      
      if(is.null(input$radioMaxGroup)){
        radio <- "2"
      } else {
        radio <- input$radioMaxGroup
      }

      chemGroupBP_group <- chemGroupBP_group()
      
      chemGroupBP_group <- chemGroupBP_group %>%
        group_by(site, date, cat) %>%
        summarise(sumEAR = sum(EAR)) %>%
        mutate(grouping=factor(cat, levels=unique(cat))) %>%  
        filter(site %in% siteToFind) 
      
      uniqueChoices <- as.character(unique(chemGroupBP_group$cat))
      
      siteLimits <- stationINFO %>%
        filter(shortName %in% unique(chemGroupBP_group$site)) 

      if(length(siteToFind) > 1){

        chemGroupBPAllSite <- chemGroupBP_group %>%
          group_by(site, cat) %>%
          summarise(maxEAR=max(sumEAR))
        
          sToxWS <- ggplot(chemGroupBPAllSite, aes(x=site, y=maxEAR, fill = cat)) +
            geom_bar(stat="identity") +
            theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.25, 
                                             colour=siteLimits$lakeColor)) +
            scale_x_discrete(limits=siteLimits$Station.shortname) +
            xlab("") +
            scale_fill_discrete("")  
          
          if(radio == "2"){
            sToxWS <- sToxWS + guides(fill=FALSE) 
          }
        
      } else {
        
        chemGroupBPOneSite <- chemGroupBP_group %>%
          select(-site)  
        
        levels(chemGroupBPOneSite$cat) <- uniqueChoices
        

        sToxWS <- ggplot(chemGroupBPOneSite, aes(x=date, y=sumEAR, fill = cat)) +
          geom_bar(stat="identity") +
          theme(axis.text.x=element_blank(),
                axis.ticks=element_blank())+
          xlab("Individual Samples") + 
          scale_fill_discrete(drop=FALSE) +
          labs(fill="")
        
        if(radio == "2"){
          sToxWS <- sToxWS + guides(fill=FALSE) 
        }

      }
      
      sToxWS
    })
    
    topPlots <- reactive({
      
      chemGroup <- chemGroup()

      if(is.null(input$group)){
        group <- unique(endPointInfo[,20])[3]
      } else {
        group <- input$group
      }
      
      if(is.null(input$sites) | input$sites == "All"){
        siteToFind <- summaryFile$site
      } else if (input$sites == "Potential 2016"){
        siteToFind <- df2016$Site
      } else {
        siteToFind <- input$sites
      }
      
      if(is.null(input$radio)){
        radio <- "2"
      } else {
        radio <- input$radio
      }
      
      chemGroupBP <- chemGroupBP()
      
      chemGroupBP_group <- chemGroupBP %>%
        filter(choices %in% group) %>%
        group_by(site, date, cat) %>%
        summarise(sumEAR = sum(EAR)) %>%
        mutate(grouping=factor(cat, levels=unique(cat))) %>%
        filter(site %in% siteToFind)
      
      uniqueChoices <- as.character(unique(chemGroupBP_group$cat))
      
      siteLimits <- stationINFO %>%
        filter(shortName %in% unique(chemGroupBP_group$site))
      
      if(length(siteToFind) > 1){
        
        chemGroupBPAllSite <- chemGroupBP_group %>%
          group_by(site, cat) %>%
          summarise(maxEAR=max(sumEAR)) %>%
          filter(site %in% siteToFind)
        
        sToxWS <- ggplot(chemGroupBPAllSite, aes(x=site, y=maxEAR, fill = cat)) +
          geom_bar(stat="identity") +
          theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.25, 
                                           colour=siteLimits$lakeColor)) +
          scale_x_discrete(limits=siteLimits$Station.shortname) +
          xlab("") +
          scale_fill_discrete("")             
          if(radio == "1"){
            sToxWS <- sToxWS + guides(fill=FALSE) 
          }
        
      } else {
        
        chemGroupBPOneSite <- chemGroupBP_group %>%
          select(-site)  
        
        levels(chemGroupBPOneSite$cat) <- uniqueChoices
        
        
        sToxWS <- ggplot(chemGroupBPOneSite, aes(x=date, y=sumEAR, fill = cat)) +
          geom_bar(stat="identity") +
          theme(axis.text.x=element_blank(),
                axis.ticks=element_blank())+
          xlab("Individual Samples") + 
          scale_fill_discrete(drop=FALSE) +
          labs(fill="")
        
          if(radio == "1"){
            sToxWS <- sToxWS + guides(fill=FALSE) 
          }       
      }
      
      sToxWS
    })
    
    bottomPlotsGroups <- reactive({
      
      if(is.null(input$sites) | input$sites == "All"){
        siteToFind <- summaryFile$site
      } else if (input$sites == "Potential 2016"){
        siteToFind <- df2016$Site
      } else {
        siteToFind <- input$sites
      }
      
      chemGroupBP_group <- chemGroupBP_group()
      
      boxData <- chemGroupBP_group %>%
        filter(site %in% siteToFind) %>%
        filter(EAR > 0)
      
      sToxWS <- ggplot(boxData) 
      
      if(!is.null(input$data) && input$data != "Passive Samples"){
        sToxWS <- sToxWS + geom_boxplot(aes(x=cat, y=EAR, fill=cat))
      } else {
        sToxWS <- sToxWS + geom_point(aes(x=cat, y=EAR, color=cat))
      }
      
      sToxWS <- sToxWS + theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.25), 
                                              legend.position = "none")+
                       xlab("") +
                       scale_y_log10("EAR") +
                       geom_hline(yintercept=0.1)
      
      sToxWS
    })

#############################################################    
    
    output$mymap <- leaflet::renderLeaflet({
      
      leaflet() %>%
        addProviderTiles("Esri.WorldPhysical") %>%
        setView(lng = -83.5, lat = 44.5, zoom=5) 
      
    })
    
    observe({
      
      if(is.null(input$sites) | input$sites == "All"){
        siteToFind <- summaryFile$site
      } else if (input$sites == "Potential 2016"){
        siteToFind <- df2016$Site
      } else {
        siteToFind <- input$sites
      }
      
      endpointSummary <- endpointSummary()
      
      sumStat <- endpointSummary %>%
        group_by(site) %>%
        summarise(meanEAR = mean(sumEAR),
                  maxEAR = max(sumEAR),
                  freq = sum(nHits > 0)/n(),
                  sumHits = sum(nHits),
                  nSamples = n()) %>%
        data.frame()%>%
        mutate(site = siteKey[site]) %>%
        filter(site %in% siteToFind)

      if(nrow(sumStat) == 0){
        sumStat <- summaryFile
        sumStat$sumHits <- sumStat$nChem
      }
      
      mapData <- right_join(stationINFO[,c("shortName", "Station.Name", "dec.lat.va","dec.long.va")], sumStat, by=c("shortName"="site"))
      mapData <- mapData[!is.na(mapData$dec.lat.va),]
      # mapData <- mapData[!is.na(mapData$nChem),]
      
      col_types <- c("darkblue","dodgerblue","green","yellow","orange","red","brown")
      leg_vals <- unique(as.numeric(quantile(mapData$maxEAR, probs=c(0,0.01,0.1,0.25,0.5,0.75,0.9,.99,1), na.rm=TRUE)))
      
      cols <- colorNumeric(col_types, domain = leg_vals)
      rad <- 1.5*seq(5000,20000, 1000)
      mapData$sizes <- rad[as.numeric(cut(mapData$nSamples, breaks=16))]
      pal = colorBin(col_types, mapData$maxEAR, bins = leg_vals)
      
      leafletProxy("mymap", data=mapData) %>%
        clearShapes() %>%
        clearControls() %>%
        addCircles(lat=~dec.lat.va, lng=~dec.long.va, 
                   popup=paste0('<b>',mapData$Station.Name,"</b><br/><table>",
                                "<tr><td>maxEAR: </td><td>",sprintf("%1.1f",mapData$maxEAR),'</td></tr>',
                                "<tr><td>nSamples: </td><td>",mapData$nSamples,'</td></tr></table>') ,
                   fillColor = ~pal(maxEAR), 
                   weight = 1,
                   color = "black",
                   fillOpacity = 0.8, radius = ~sizes, opacity = 0.8) %>%
        addLegend(
          position = 'bottomleft',
          pal=pal,
          values=~maxEAR,
          opacity = 0.8,
          title = 'Maximum EAR')
      
    })
    
############################################################# 
    
    output$groupControl <- renderUI({

#       if(is.null(input$groupCol)){
#         groupCol <- names(endPointInfo)[20]
#       } else {
#         groupCol <- input$groupCol
#       }
#       
#       statCol <- statsOfColumn()
#       
#       freqCol <- grep("freq",names(statCol))
#       maxEARS <- grep("maxEAR",names(statCol))
#       
#       statCol <- statCol[,c(1,c(maxEARS,freqCol)[order(c(maxEARS,freqCol))])]
#       
#       maxEARS <- grep("maxEAR",names(statCol))
#       
#       MaxEARSordered <- order(apply(statCol[,maxEARS], 2, max),decreasing = TRUE)
#       
#       statCol <- statCol[,c(1,interl(maxEARS[MaxEARSordered],(maxEARS[MaxEARSordered]-1)))]
#       
#       freqCol <- grep("freq",names(statCol))
#       maxEARS <- grep("maxEAR",names(statCol))
#       
#       namesToUse <- gsub("maxEAR","",names(statCol)[-1])
#       namesToUse <- gsub("freq","",namesToUse)
#       namesToUse <- unique(namesToUse)
#       namesToUse <- gsub("^\\s+|\\s+$", "", namesToUse)
# 
#       nEndPointsInChoice <- as.character(table(endPointInfo[,input$groupCol])[namesToUse])
#       dropDownHeader <- paste0(namesToUse," (",nEndPointsInChoice,")")
      
      ChoicesInGroup <- names(table(endPointInfo[,input$groupCol]))
      nEndPointsInChoice <- as.character(table(endPointInfo[,input$groupCol]))
      dropDownHeader <- paste0(ChoicesInGroup," (",nEndPointsInChoice,")")
      
      selectInput("group", label = "Group in annotation (# End Points)",
                  choices = setNames(ChoicesInGroup,dropDownHeader),
                  multiple = FALSE)

    })
    
    output$TableHeader <- renderUI({
      HTML(paste("<h4>", input$group,"-",input$data,": ",input$sites, "</h4>"))
    })
    
    output$BoxHeader <- renderUI({
      HTML(paste("<h4>", input$group,"-",input$data,": ",input$sites, "</h4>"))
    })

    output$TableHeaderColumns <- renderUI({
      HTML(paste("<h4>", input$data,": ",input$sites, "</h4>"))
    })

    output$nGroup <- renderUI({
      
      if(is.null(input$radioMaxGroup)){
        radio <- 1
      } else {
        radio <- input$radioMaxGroup
      }
      if(radio == "2"){
          word <- "chemicals"
      } else if (radio == "3") {
        word <- "classes"
      } else {
        word <- "choices"
      }
      
      textUI <- paste("<h5>max = Maximum of number of",word,"with hits per sample</h5>",
                      "<h5>mean = Mean number of",word,"with hits per sample</h5>",
                      "<h5>count = Number of",word,"with hits per sample</h5>")
      
      HTML(textUI)
      
    })

})