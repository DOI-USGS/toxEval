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
# wData <- wData
# pCodeInfo <- pCodeInfo

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

choicesPerGroup <- apply(endPointInfo[,-3], 2, function(x) length(unique(x)))
groupChoices <- paste0(names(choicesPerGroup)," (",choicesPerGroup,")")


shinyServer(function(input, output) {
  
    
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
      
      chemicalSummary
      
    })
  
    chemicalSummaryFiltered <- reactive({
      
      if(is.null(input$groupCol)){
        groupCol <- names(endPointInfo)[20]
      } else {
        groupCol <- input$groupCol
      }
      
      endPointInfoSub <- select_(endPointInfo, "assay_component_endpoint_name", groupCol) %>%
        distinct()
      
      chemicalSummaryFiltered <- chemicalSummary() %>%
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
    
    chemSiteSumm <- reactive({
      
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
      
      if(is.null(input$sites) | input$sites == "All"){
        siteToFind <- summaryFile$site
      } else {
        siteToFind <- input$sites
      }
      
      chemSiteSumm <- chemicalSummaryFiltered %>%
        filter_(paste0(groupCol," == '", group, "'")) %>%
        mutate(site = siteKey[site]) %>%
        filter_(paste0("site == '", siteToFind, "'")) %>%
        group_by(chnm, date) %>%
        summarise(sumEAR = sum(EAR),
                  nHits = sum(hits)) 
      
      chemSiteSumm
    })
    
    groupSiteSumm <- reactive({
      
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
      
      if(is.null(input$sites) | input$sites == "All"){
        siteToFind <- summaryFile$site
      } else {
        siteToFind <- input$sites
      }
      
      groupSiteSumm <- chemicalSummaryFiltered %>%
        filter_(paste0(groupCol," == '", group, "'")) %>%
        mutate(site = siteKey[site]) %>%
        filter_(paste0("site == '", siteToFind, "'")) %>%
        group_by(chnm, date) %>%
        summarise(sumEAR = sum(EAR),
                  nHits = sum(hits)) 
      
      groupSiteSumm
    })
    
    endpointSummary <- reactive({
      
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
      
      endpointSummary <- chemicalSummaryFiltered() %>%
        filter_(paste0(groupCol," == '", group, "'")) %>%
        group_by(site, date) %>%
        summarise(sumEAR = sum(EAR),
                  nHits = sum(hits)) 
    })
    
    statsOfSum <- reactive({
      endpointSummary() %>%
        group_by(site) %>%
        summarise(meanEAR = mean(sumEAR),
                  maxEAR = max(sumEAR),
                  sumHits = sum(nHits),
                  nSamples = n()) %>%
        data.frame()%>%
        mutate(site = siteKey[site]) 
      
    })
    
    statsOfChem <- reactive({
      chemSiteSumm() %>%
        group_by(chnm) %>%
        summarise(meanEAR = mean(sumEAR),
                  maxEAR = max(sumEAR),
                  sumHits = sum(nHits)) %>%
        data.frame()
      
    })
    
    statsOfColumn <- reactive({
      if(is.null(input$groupCol)){
        groupCol <- names(endPointInfo)[20]
      } else {
        groupCol <- input$groupCol
      }
      
      if(input$sites == "All"){
        
        statsOfColumn <- chemicalSummary() %>%
          rename(assay_component_endpoint_name=endPoint) %>%
          filter(assay_component_endpoint_name %in% endPointInfo$assay_component_endpoint_name ) %>%
          data.table()%>%
          left_join(data.table(endPointInfo[,c("assay_component_endpoint_name", groupCol)]), by = "assay_component_endpoint_name") %>%
          data.frame()%>%
          rename(endPoint=assay_component_endpoint_name)%>%
          select_("hits","EAR","endPoint","site","date","choices"=groupCol) %>%
          group_by(site, date,choices) %>%
          summarise(sumEAR = sum(EAR),
                    nHits = sum(hits)) %>%
          group_by(site,choices) %>%
          summarise(maxEAR = max(sumEAR),
                    sumHits = sum(nHits)) %>%
          data.frame()%>%
          mutate(site = siteKey[site]) %>%
          gather(calc, value, -site, -choices) %>%
          unite(choice_calc, choices, calc, sep=" ") %>%
          spread(choice_calc, value)
      } else {
        
        if(is.null(input$sites)){
          siteToFind <- summaryFile$site
        } else {
          siteToFind <- input$sites
        }
        
        statsOfColumn <- chemicalSummary() %>%
          rename(assay_component_endpoint_name=endPoint) %>%
          filter(assay_component_endpoint_name %in% endPointInfo$assay_component_endpoint_name ) %>%
          data.table()%>%
          left_join(data.table(endPointInfo[,c("assay_component_endpoint_name", groupCol)]), by = "assay_component_endpoint_name") %>%
          data.frame() %>%
          mutate(site = siteKey[site]) %>%
          filter_(paste0("site == '", siteToFind, "'")) %>%
          select_("hits","EAR","chnm","class","date","choices"=groupCol) %>%
          group_by(chnm, date,choices) %>%
          summarise(sumEAR = sum(EAR),
                    nHits = sum(hits)) %>%
          group_by(chnm,choices) %>%
          summarise(maxEAR = max(sumEAR),
                    sumHits = sum(nHits)) %>%
          data.frame() %>%
          gather(calc, value, -chnm, -choices) %>%
          unite(choice_calc, choices, calc, sep=" ") %>%
          spread(choice_calc, value)
      }
      
    })
    
    chemGroup <- reactive({
      
      if(is.null(input$groupCol)){
        groupCol <- names(endPointInfo)[20]
      } else {
        groupCol <- input$groupCol
      }
      
      chemGroup <- chemicalSummary() %>%
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
      
      statsOfGroup <- chemGroup()
      
      if(input$sites == "All"){
        statsOfGroup <-  statsOfGroup %>%
          group_by(site, date,choices) %>%
          summarise(nChemWithHits = length(unique(chnm[EAR > 0.1]))) %>%
          group_by(site,choices) %>%
          summarise(maxChem = max(nChemWithHits),
                    meanChem = mean(nChemWithHits),
                    sumChemWithHits = sum(nChemWithHits)) %>%
          data.frame()%>%
          gather(calc, value, -site, -choices) %>%
          unite(choice_calc, choices, calc, sep=" ") %>%
          spread(choice_calc, value)
      } else {
        
        if(is.null(input$sites)){
          siteToFind <- summaryFile$site
        } else {
          siteToFind <- input$sites
        }
        
        statsOfGroup <-  statsOfGroup %>%
          filter_(paste0("site == '", siteToFind, "'")) %>%
          group_by(site, date,choices) %>%
          summarise(nChemWithHits = length(unique(chnm[EAR > 0.1]))) %>%
          group_by(site,choices) %>%
          summarise(maxChem = max(nChemWithHits),
                    meanChem = mean(nChemWithHits),
                    sumChemWithHits = sum(nChemWithHits)) %>%
          data.frame() %>%
          select(-site)
      }
    })
    
    output$groupControl <- renderUI({
      
      ChoicesInGroup <- names(table(endPointInfo[,input$groupCol]))
      nEndPointsInChoice <- as.character(table(endPointInfo[,input$groupCol]))
      dropDownHeader <- paste0(ChoicesInGroup," (",nEndPointsInChoice,")")
      
      selectInput("group", label = "Group in annotation (# End Points)",
                  choices = setNames(ChoicesInGroup,dropDownHeader),
                  selected = unique(endPointInfo[,20])[3],
                  multiple = FALSE)
    })
    
    output$TableHeader <- renderUI({
      HTML(paste("<br/><br/><h3>Table of summations summaries:",input$group,"-",input$data,"</h3>"))
    })
    
    output$BoxHeader <- renderUI({
      HTML(paste("<br/><h3>Boxplot summaries:",input$group,"-",input$data,"</h3>"))
    })
    
    output$TableHeaderColumns <- renderUI({
      HTML(paste("<br/><h3>Table of summations summaries:",input$groupCol,"-",input$data,"</h3>"))
    })
    
    output$TableHeaderColumns2 <- renderUI({
      HTML(paste("<br/><h3>Table of chemical summaries:",input$groupCol,"-",input$data,"</h3>"))
    })
    
    output$BoxHeaderColumns <- renderUI({
      HTML(paste("<br/><h3>Boxplot summaries:",input$groupCol,"-",input$data,"</h3>"))
    })
    
    output$table <- DT::renderDataTable({
      
      if(input$sites == "All"){
        statsOfSumDF <- DT::datatable(statsOfSum(), 
                                      rownames = FALSE,
                                      colnames = c('Maximum EAR' = 3, 'Sum of Hits' = 4),
                                      filter = 'top',
                                      options = list(pageLength = 10, 
                                                     order=list(list(2,'desc')))) %>%
          formatRound(c("Maximum EAR","meanEAR"), 2) %>%
          formatStyle("Maximum EAR", 
                      background = styleColorBar(statsOfSum()[,3], 'steelblue'),
                      backgroundSize = '100% 90%',
                      backgroundRepeat = 'no-repeat',
                      backgroundPosition = 'center'
          ) %>%
          formatStyle("Sum of Hits", 
                      background = styleColorBar(statsOfSum()[,4], 'wheat'),
                      backgroundSize = '100% 90%',
                      backgroundRepeat = 'no-repeat',
                      backgroundPosition = 'center'
          )
      } else {
        statsOfSumDF <- DT::datatable(statsOfChem(), 
                                      rownames = FALSE,
                                      colnames = c('Chemical Name' = 1, 'Mean EAR' = 2,
                                                   'Maximum EAR' = 3, 'Sum of Hits' = 4),
                                      filter = 'top',
                                      options = list(pageLength = 10, 
                                                     order=list(list(2,'desc')))) %>%
          formatRound(c("Maximum EAR","Mean EAR"), 2) %>%
          formatStyle("Maximum EAR", 
                      background = styleColorBar(statsOfChem()[,3], 'steelblue'),
                      backgroundSize = '100% 90%',
                      backgroundRepeat = 'no-repeat',
                      backgroundPosition = 'center'
          ) %>%
          formatStyle("Sum of Hits", 
                      background = styleColorBar(statsOfChem()[,4], 'wheat'),
                      backgroundSize = '100% 90%',
                      backgroundRepeat = 'no-repeat',
                      backgroundPosition = 'center'
          )
      }
    })
    
    output$tableGroupSumm <- DT::renderDataTable({
      
      statsOfGroup <- statsOfGroup() 
      
      if(input$sites == "All"){
        
        maxChem <- grep("maxChem",names(statsOfGroup))
        
        MaxChemordered <- order(apply(statsOfGroup[,maxChem], 2, max),decreasing = TRUE)
        
        
        interl3 <- function (a,b,cv) {
          n <- min(length(a),length(b),length(cv))
          p1 <- as.vector(rbind(a[1:n],b[1:n],cv[1:n]))
          p2 <- c(a[-(1:n)],b[-(1:n)],cv[-(1:n)])
          c(p1,p2)
        }
        
        if(length(maxChem) > 9){
          statsOfGroup <- statsOfGroup[,c(1,interl3(maxChem[MaxChemordered[1:9]],(maxChem[MaxChemordered[1:9]]+1),(maxChem[MaxChemordered[1:9]]+2)))]
          maxChem <- maxChem[1:9]
        } else {
          statsOfGroup <- statsOfGroup[,c(1,interl3(maxChem[MaxChemordered],(maxChem[MaxChemordered]+1),(maxChem[MaxChemordered]+2)))]
        }
        
        colors <- brewer.pal(length(maxChem),"Blues") #"RdYlBu"
        
        tableGroup <- DT::datatable(statsOfGroup, 
                                    rownames = FALSE,
                                    filter = 'top',
                                    options = list(pageLength = 10, 
                                                   order=list(list(1,'desc'))))
        tableGroup <- formatRound(tableGroup, names(statsOfGroup)[grep("meanChem",names(statsOfGroup))], 2)
        
        
        for(i in 1:length(maxChem)){
          tableGroup <- formatStyle(tableGroup, 
                                    names(statsOfGroup)[maxChem[i]], 
                                    backgroundColor = colors[i])
          tableGroup <- formatStyle(tableGroup, 
                                    names(statsOfGroup)[maxChem[i]+1], 
                                    backgroundColor = colors[i])
          tableGroup <- formatStyle(tableGroup, 
                                    names(statsOfGroup)[maxChem[i]+2], 
                                    backgroundColor = colors[i])
          
          tableGroup <- formatStyle(tableGroup, names(statsOfGroup)[maxChem[i]], 
                                    background = styleColorBar(statsOfGroup[,maxChem[i]], 'goldenrod'),
                                    backgroundSize = '100% 90%',
                                    backgroundRepeat = 'no-repeat',
                                    backgroundPosition = 'center' ) 
          
          tableGroup <- formatStyle(tableGroup, names(statsOfGroup)[maxChem[i]+1], 
                                    background = styleColorBar(statsOfGroup[,maxChem[i]+1], 'wheat'),
                                    backgroundSize = '100% 90%',
                                    backgroundRepeat = 'no-repeat',
                                    backgroundPosition = 'center') 
          
          tableGroup <- formatStyle(tableGroup, names(statsOfGroup)[maxChem[i]+2], 
                                    background = styleColorBar(statsOfGroup[,maxChem[i]+2], 'seashell'),
                                    backgroundSize = '100% 90%',
                                    backgroundRepeat = 'no-repeat',
                                    backgroundPosition = 'center') 
        }
      } else {
        tableGroup <- DT::datatable(statsOfGroup, 
                                    rownames = FALSE,
                                    filter = 'top',
                                    colnames = c('Annotation' = 1, 'Maximum Number of Chemicals per Sample with Hits' = 2,
                                                 'Mean Number of Chemicals per Sample with Hits' = 3, 'Total Number of Chemicals with Hits' = 4),
                                    options = list(pageLength = 10, 
                                                   order=list(list(1,'desc'))))
        tableGroup <- formatRound(tableGroup, 'Mean Number of Chemicals per Sample with Hits', 2)
        tableGroup <- formatStyle(tableGroup, 'Maximum Number of Chemicals per Sample with Hits', 
                                  background = styleColorBar(statsOfGroup[,grep("maxChem",names(statsOfGroup))], 'goldenrod'),
                                  backgroundSize = '100% 90%',
                                  backgroundRepeat = 'no-repeat',
                                  backgroundPosition = 'center' ) 
        
        tableGroup <- formatStyle(tableGroup, 'Mean Number of Chemicals per Sample with Hits', 
                                  background = styleColorBar(statsOfGroup[,grep("meanChem",names(statsOfGroup))], 'wheat'),
                                  backgroundSize = '100% 90%',
                                  backgroundRepeat = 'no-repeat',
                                  backgroundPosition = 'center') 
        
        tableGroup <- formatStyle(tableGroup, 'Total Number of Chemicals with Hits', 
                                  background = styleColorBar(statsOfGroup[,grep("sumChemWithHits",names(statsOfGroup))], 'seashell'),
                                  backgroundSize = '100% 90%',
                                  backgroundRepeat = 'no-repeat',
                                  backgroundPosition = 'center')
      }
      
      tableGroup
      
    })
    
    output$tableSumm <- DT::renderDataTable({
      statCol <- statsOfColumn()
      
      maxEARS <- grep("maxEAR",names(statCol))
      
      MaxEARSordered <- order(apply(statCol[,maxEARS], 2, max),decreasing = TRUE)
      
      interl <- function (a,b) {
        n <- min(length(a),length(b))
        p1 <- as.vector(rbind(a[1:n],b[1:n]))
        p2 <- c(a[-(1:n)],b[-(1:n)])
        c(p1,p2)
      }
      
      if(length(maxEARS) > 9){
        statCol <- statCol[,c(1,interl(maxEARS[MaxEARSordered[1:9]],(maxEARS[MaxEARSordered[1:9]]+1)))]
        maxEARS <- maxEARS[1:9]
      } else {
        statCol <- statCol[,c(1,interl(maxEARS[MaxEARSordered],(maxEARS[MaxEARSordered]+1)))]
      }
      
      colors <- brewer.pal(length(maxEARS),"Blues") #"RdYlBu"
      tableSumm <- DT::datatable(statCol, 
                                 rownames = FALSE,
                                 options = list(order=list(list(1,'desc'))))
      
      tableSumm <- formatRound(tableSumm, names(statCol)[maxEARS], 2) 
      
      for(i in 1:length(maxEARS)){
        tableSumm <- formatStyle(tableSumm, 
                                 names(statCol)[maxEARS[i]], 
                                 backgroundColor = colors[i])
        tableSumm <- formatStyle(tableSumm, 
                                 names(statCol)[maxEARS[i]+1], 
                                 backgroundColor = colors[i])
        
        tableSumm <- formatStyle(tableSumm, names(statCol)[maxEARS[i]], 
                                 background = styleColorBar(statCol[,maxEARS[i]], 'goldenrod'),
                                 backgroundSize = '100% 90%',
                                 backgroundRepeat = 'no-repeat',
                                 backgroundPosition = 'center' ) 
        tableSumm <- formatStyle(tableSumm, names(statCol)[maxEARS[i]+1], 
                                 background = styleColorBar(statCol[,maxEARS[i]+1], 'wheat'),
                                 backgroundSize = '100% 90%',
                                 backgroundRepeat = 'no-repeat',
                                 backgroundPosition = 'center') 
        
      }
      
      tableSumm
      
    })
    
    output$stackBar <- renderPlot({ 
      print(topPlots())
    })
    
    output$graph <- renderPlot({ 
      print(bottomPlots())
    })
    
    topPlots <- reactive({
      chemGroup <- chemGroup()
      
      if(is.null(input$sites)){
        siteToFind <- summaryFile$site
      } else {
        siteToFind <- input$sites
      }
      
      if(is.null(input$group)){
        group <- unique(endPointInfo[,20])[3]
      } else {
        group <- input$group
      }
      
      if(nrow(chemGroup) == 0){
        chemGroup$site <- stationINFO$fullSiteID
      }
      
      siteLimits <- stationINFO %>%
        filter(shortName %in% unique(chemGroup$site))
      
      chemGroupBP <- chemGroup %>%
        filter(EAR >= 0.1) %>%
        select(EAR, chnm, site, choices, date) %>%
        data.frame() %>%
        filter_(paste0("choices == '", group, "'")) %>%
        select(-choices) %>%
        mutate(chnm=factor(chnm, levels=unique(chnm))) %>%
        mutate(date = factor(as.numeric(date) - min(as.numeric(date))))%>%
        arrange(chnm)
      
      uniqueChms <- as.character(unique(chemGroupBP$chnm))
      
      if(siteToFind == "All"){
        
        sToxWS <- ggplot(chemGroupBP, aes(x=site, y=EAR, fill = chnm)) +
          geom_bar(stat="identity") +
          theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.25, 
                                           colour=siteLimits$lakeColor)) +
          scale_x_discrete(limits=siteLimits$Station.shortname) +
          xlab("")             
        
      } else {
        
        chemGroupBPOneSite <- chemGroupBP %>%
          filter_(paste0("site == '", siteToFind, "'")) %>%
          select(-site)  
        
        levels(chemGroupBPOneSite$chnm) <- uniqueChms
        
        sToxWS <- ggplot(chemGroupBPOneSite, aes(x=date, y=EAR, fill = chnm)) +
          geom_bar(stat="identity") +
          theme(axis.text.x=element_blank(),
                axis.ticks=element_blank())+
          xlab("Individual Samples") + 
          scale_fill_discrete(drop=FALSE)
        
      }
      
      sToxWS
    })
    
    bottomPlots <- reactive({
      
      g <- ggplot_build(topPlots())
      fillColors <- g$data[[1]]["fill"][[1]]
      chnms <- g$plot$data$chnm
      
      dfNames <- data.frame(chnm = names(table(chnms)[table(chnms) != 0]),
                            freq = as.numeric(table(chnms)[table(chnms) != 0]),stringsAsFactors = FALSE)
      dfCol <- data.frame(colors = names(table(fillColors)), 
                 freq = as.numeric(table(fillColors)),stringsAsFactors = FALSE)
      
      colorKey <- left_join(dfNames, dfCol, by="freq")

      if(input$sites == "All"){
        
        endpointSummary <- endpointSummary()
        
        if(nrow(endpointSummary) == 0){
          endpointSummary$site <- stationINFO$fullSiteID
        }
        
        siteLimits <- stationINFO %>%
          filter(fullSiteID %in% unique(endpointSummary$site))
        
        endPointSummBP <- endpointSummary %>%
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
        chemSiteSumm <- chemSiteSumm() %>%
          mutate(sumEARnoZero = sumEAR) 
        
        ndLevel <- 0.1*min(chemSiteSumm$sumEARnoZero[chemSiteSumm$sumEARnoZero != 0])
        
        if(is.finite(ndLevel)){
          chemSiteSumm$sumEARnoZero[chemSiteSumm$sumEARnoZero == 0] <- ndLevel
    
          noFill <- data.frame(chnm=unique(chemSiteSumm$chnm)[!(unique(chemSiteSumm$chnm) %in% colorKey$chnm)],
                               colors=rep("#FFFFFF",length(unique(chemSiteSumm$chnm))-nrow(colorKey)),
                               stringsAsFactors = FALSE)
          fillTotal <- rbind(colorKey[,c('colors','chnm')], noFill)
          
          sToxWS <- ggplot(chemSiteSumm) 
          
          if(!is.null(input$data) && input$data != "Passive Samples"){
            sToxWS <- sToxWS + geom_boxplot(aes(x=chnm, y=sumEAR, 
                                                fill=chnm)) +
              scale_fill_manual(values=setNames(fillTotal$colors,fillTotal$chnm))
          } else {
            sToxWS <- sToxWS + geom_point(aes(x=chnm, y=sumEAR)) 
          }
          sToxWS <- sToxWS + theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.25), 
                                   legend.position = "none")+
            xlab("") +
            scale_y_log10("Summation of EAR per sample")
        }
      }
      
      sToxWS
    })
    
    output$mymap <- leaflet::renderLeaflet({
      
      leaflet() %>%
        addProviderTiles("Esri.WorldPhysical") %>%
        setView(lng = -83.5, lat = 44.5, zoom=5) 
      
    })
    
    observe({
      
      sumStat <- statsOfSum()
      
      if(nrow(sumStat) == 0){
        sumStat <- summaryFile
        sumStat$sumHits <- sumStat$nChem
      }
      
      mapData <- right_join(stationINFO, sumStat, by=c("shortName"="site"))
      mapData <- mapData[!is.na(mapData$dec.lat.va),]
      # mapData <- mapData[!is.na(mapData$nChem),]
      
      col_types <- c("darkblue","dodgerblue","green","yellow","orange","red","brown")
      leg_vals <- unique(as.numeric(quantile(mapData$maxEAR, probs=c(0,0.01,0.1,0.25,0.5,0.75,0.9,.99,1), na.rm=TRUE)))
      
      cols <- colorNumeric(col_types, domain = leg_vals)
      rad <- 1.5*seq(5000,25000, 500)
      mapData$sizes <- rad[as.numeric(
        cut(mapData$sumHits, 
            c(-1,unique(as.numeric(
              quantile(mapData$sumHits, probs=seq(.01,.99,length=length(rad)), na.rm=TRUE)
            )),(max(mapData$sumHits,na.rm=TRUE)+1))
        ))]
      pal = colorBin(col_types, mapData$maxEAR, bins = leg_vals)
      
      if(input$sites != "All"){
        mapData <- mapData[mapData$shortName == input$sites,]
      }
      
      leafletProxy("mymap", data=mapData) %>%
        clearShapes() %>%
        clearControls() %>%
        addCircles(lat=~dec.lat.va, lng=~dec.long.va, 
                   popup=mapData$Station.Name
                   #                    paste0('<b>',mapData$Station.Name,"</b><br/><table>",
                   #                           "<tr><td>Frequency</td><td>",sprintf("%1.1f",mapData$freq),'</td></tr>',
                   #                           "<tr><td>nSample</td><td>",mapData$nSamples,'</td></tr>',
                   #                           "<tr><td>maxEAR</td><td>",sprintf("%1.1f",mapData$maxEAR),'</td></tr>',
                   #                           "<tr><td>nChem</td><td>",mapData$nChem,'</td></tr>',
                   #                           "<tr><td>nEndPoints</td><td>",mapData$nEndPoints,'</td></tr></table>')
                   ,
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

})