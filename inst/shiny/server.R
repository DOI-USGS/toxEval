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

makePlots <- function(boxData, noLegend, boxPlot, siteToFind){
  
  siteToFind <- unique(boxData$site)

  graphData <- boxData %>%
    group_by(site,date,cat) %>%
    summarise(sumEAR=sum(EAR)) %>%
    data.frame() %>%
    mutate(grouping=factor(cat, levels=unique(cat)))  %>%
    group_by(site, cat) %>%
    summarise(maxEAR=max(sumEAR),
              meanEAR=mean(sumEAR)) %>%
    gather(stat, value, -site, -cat)

  if(length(siteToFind) > 1){
    lowerPlot <- ggplot(graphData[graphData$stat == "meanEAR",])+
      scale_y_log10("Mean EAR Per Site")
  } else {
    siteData <- boxData %>%
      group_by(site,date,cat) %>%
      summarise(sumEAR=sum(EAR)) %>%
      data.frame() %>%
      mutate(grouping=factor(cat, levels=unique(cat))) %>%
      rename(value=sumEAR)
      
    lowerPlot <- ggplot(siteData)+
      scale_y_log10("Sum of EAR")
  }
   
  
  if(!boxPlot){
    lowerPlot <- lowerPlot + geom_point(aes(x=cat, y=value, color=cat, size=3))
  } else {
    lowerPlot <- lowerPlot + geom_boxplot(aes(x=cat, y=value, fill=cat))
  }
  
  lowerPlot <- lowerPlot + theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.25), 
                                 legend.position = "none") +
    xlab("") +
    geom_hline(yintercept=0.1)  

  uniqueChoices <- as.character(unique(boxData$cat))
  
  siteLimits <- stationINFO %>%
    filter(shortName %in% unique(boxData$site))
  
  if(length(siteToFind) > 1){

    upperPlot <- ggplot(graphData, aes(x=site, y=value, fill = cat)) +
      geom_bar(stat="identity") +
      facet_wrap(~stat, nrow=2, ncol=1, scales = "free_y") + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.25, 
                                       colour=siteLimits$lakeColor)) +
      scale_x_discrete(limits=siteLimits$Station.shortname) +
      xlab("") +
      ylab("EAR") +
      scale_fill_discrete("") 
      
    
    if(noLegend){
      upperPlot <- upperPlot + guides(fill=FALSE) 
    }

  } else {
    
    chemGroupBPOneSite <- boxData %>%
      select(-site)  
    
    levels(chemGroupBPOneSite$cat) <- uniqueChoices
    
    upperPlot <- ggplot(chemGroupBPOneSite, aes(x=date, y=EAR, fill = cat)) +
      geom_bar(stat="identity")+
      theme(axis.text.x=element_blank(),
            axis.ticks=element_blank())+
      xlab("Individual Samples") + 
      ylab("EAR") +
      scale_fill_discrete(drop=FALSE) +
      labs(fill="") 
    
    if(noLegend){
      upperPlot <- upperPlot + guides(fill=FALSE) 
    } 
  }
  
  return(list(upperPlot=upperPlot, lowerPlot=lowerPlot))
}

shinyServer(function(input, output) {
  
#############################################################   
    chemicalSummary <- reactive({
      
      if(is.null(input$data) | input$data == "Water Sample"){
        chemicalSummary <- chemicalSummaryWS
      } else if (input$data == "Passive Samples"){
        chemicalSummary <- chemicalSummaryPS
      } else {
        chemicalSummary <- rbind(chemicalSummaryWS, chemicalSummaryPS)
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
        rename(endPoint = assay_component_endpoint_name) %>%
        select_("hits","EAR","chnm","class","date","choices"=groupCol,"site","endPoint","endPointValue") %>%
        filter(site %in% siteToFind)
      
      if(length(siteToFind) == 1){
 
        statsOfColumn <- statsOfColumn %>%
          # mutate(site=siteKey[site]) %>%
          filter(site %in% siteToFind) 
        
        if(radio == "2"){
          statsOfColumn <- mutate(statsOfColumn, site = chnm) 
        } else if (radio == "3"){
          statsOfColumn <- mutate(statsOfColumn, site = class) 
        } else if (radio == "4"){
          statsOfColumn <- mutate(statsOfColumn, site = endPoint)
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
                  meanEAR = mean(sumEAR),
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
        select_("hits","EAR","chnm","class","site","date","choices"=groupCol,"endPoint", "endPointValue")
      
    })
    
    statsOfGroupOrdered <- reactive({
      
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
      } else if(radioMaxGroup == "4"){
        statsOfGroup <- statsOfGroup %>%
          mutate(category=endPoint) 
      } else {
        statsOfGroup <- statsOfGroup %>%
          mutate(category = choices)         
      }
      
      statsOfGroupOrdered <- statsOfGroup %>%
        group_by(site, date,category) %>%
        summarise(sumEAR = sum(EAR)) %>%
        group_by(site,category) %>%
        summarise(max = sum(sumEAR > 0.1),
                  mean = sum(mean(sumEAR) > 0.1),
                  hit = as.numeric(any(sumEAR > 0.1)) )%>%
        data.frame()

      if(length(siteToFind) > 1){
        statsOfGroupOrdered <- statsOfGroupOrdered %>%
          group_by(site) %>%
          summarise(max=sum(max > 0),
                    mean=sum(mean > 0))
      } else {
        statsOfGroupOrdered <- statsOfGroupOrdered %>%
          data.frame() %>%
          left_join(summaryFile[c("site","nSamples")], by="site") %>%
          select(-site, -hit)
      }
      
      statsOfGroupOrdered
        
    }) 
    
############################################################# 
        
    output$table <- DT::renderDataTable({

      if(is.null(input$radio)){
        radio <- "1"
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
      } else if(radio == "3"){
        chemicalSummaryFiltered <- chemicalSummaryFiltered %>%
          mutate(grouping=endPoint)  
      } else {
        chemicalSummaryFiltered <- chemicalSummaryFiltered %>%
          mutate(grouping=chnm)        
      }
      
      statTable <- chemicalSummaryFiltered %>%
        filter_(paste0(groupCol," == '", group, "'")) %>%
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

      if(length(siteToFind) > 1){
        
        if (input$data == "Water Sample"){
          
          summaryFile <- filter(summaryFile, site %in% siteToFind)
          
          statsOfGroupOrdered <- left_join(summaryFile[,c("site","nSamples")], statsOfGroupOrdered, by="site")
        }
        
        meanChem <- grep("mean",names(statsOfGroupOrdered))
        maxChem <- grep("max",names(statsOfGroupOrdered))

        tableGroup <- DT::datatable(statsOfGroupOrdered, 
                                      rownames = FALSE,
                                      filter = 'top',
                                      options = list(pageLength = nrow(statsOfGroupOrdered), 
                                                     order=list(list(2,'desc'))))

        tableGroup <- formatStyle(tableGroup, names(statsOfGroupOrdered)[maxChem], 
                                  background = styleColorBar(range(statsOfGroupOrdered[,maxChem],na.rm=TRUE), 'goldenrod'),
                                  backgroundSize = '100% 90%',
                                  backgroundRepeat = 'no-repeat',
                                  backgroundPosition = 'center' ) 
        
        tableGroup <- formatStyle(tableGroup, names(statsOfGroupOrdered)[meanChem], 
                                  background = styleColorBar(range(statsOfGroupOrdered[,meanChem],na.rm=TRUE), 'wheat'),
                                  backgroundSize = '100% 90%',
                                  backgroundRepeat = 'no-repeat',
                                  backgroundPosition = 'center') 

        } else {
          
          tableGroup <- DT::datatable(statsOfGroupOrdered[,c("category","max","nSamples")], 
                                      rownames = FALSE,
                                      colnames = c('hits' = 2),
                                      filter = 'top',
                                      options = list(pageLength = nrow(statsOfGroupOrdered), 
                                                     order=list(list(1,'desc'))))

          tableGroup <- formatStyle(tableGroup, "hits", 
                                    background = styleColorBar(range(statsOfGroupOrdered[,"max"],na.rm=TRUE), 'goldenrod'),
                                    backgroundSize = '100% 90%',
                                    backgroundRepeat = 'no-repeat',
                                    backgroundPosition = 'center' ) 

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
      print(groupPlots()$upperPlot)
    })

    output$stackBar <- renderPlot({ 
      print(plots()$upperPlot)
    })

    output$graph <- renderPlot({ 
      print(plots()$lowerPlot)
    })
    
    output$graphGroup <- renderPlot({ 
      print(groupPlots()$lowerPlot)
    })
    
    boxData2 <- reactive({
      
      if(is.null(input$radio)){
        noLegend <- TRUE
        radio <- "1"
      } else {
        noLegend <- input$radio %in% c("1","3")
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

      chemGroup <- chemGroup()
      
      chemGroupBP_group <- chemGroup %>%
        mutate(date = factor(as.numeric(date) - min(as.numeric(date))))
      
      if(radio == "1"){
        chemGroupBP_group <- mutate(chemGroupBP_group, cat=chnm)
      } else if (radio == "2"){
        chemGroupBP_group <- mutate(chemGroupBP_group, cat=class)
      } else if (radio == "3"){
        chemGroupBP_group <- mutate(chemGroupBP_group, cat=endPoint) %>%
          filter(EAR > 0.1) 
      }
      
      boxData <- chemGroupBP_group %>%
        filter(site %in% siteToFind) %>%
        filter(choices %in% group) 
        # filter(EAR > 0)
      
      return(boxData)
      
    })
    
    plots <- reactive({
      
      boxData <- boxData2()
      
      if(is.null(input$sites) | input$sites == "All"){
        siteToFind <- summaryFile$site
      } else if (input$sites == "Potential 2016"){
        siteToFind <- df2016$Site
      } else {
        siteToFind <- input$sites
      }

      if(is.null(input$radio)){
        noLegend <- TRUE
      } else {
        noLegend <- input$radio %in% c("1","3")
      }
      
      if(is.null(input$data)){
        boxPlot <- TRUE
      } else {
        boxPlot <- !(input$data == "Passive Samples" & length(siteToFind) == 1) & !(input$radio == "3")
      }
      
      return(makePlots(boxData, noLegend, boxPlot, siteToFind))
      
    })

    boxData <- reactive({
      
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
      
      chemGroupBP_group <- chemGroup %>%
        mutate(date = factor(as.numeric(date) - min(as.numeric(date))))
      
      if(radioMaxGroup == "1"){
        chemGroupBP_group <- mutate(chemGroupBP_group, cat=choices)
      } else if (radioMaxGroup == "2"){
        chemGroupBP_group <- mutate(chemGroupBP_group, cat=chnm)
      } else if (radioMaxGroup == "4"){
        chemGroupBP_group <- mutate(chemGroupBP_group, cat=endPoint) %>%
          filter(EAR > 0.1) 
      } else {
        chemGroupBP_group <- mutate(chemGroupBP_group, cat=class)
      }
      
      boxData <- chemGroupBP_group %>%
        filter(site %in% siteToFind)
      
      return(boxData)
      
    })
    
    groupPlots <- reactive({
      
      boxData <- boxData()
      
      if(is.null(input$sites) | input$sites == "All"){
        siteToFind <- summaryFile$site
      } else if (input$sites == "Potential 2016"){
        siteToFind <- df2016$Site
      } else {
        siteToFind <- input$sites
      }

      noLegend <- FALSE
      radioMaxGroup <- "1"
      boxGraph <- TRUE
      
      if(!is.null(input$radioMaxGroup)){
        noLegend <- input$radioMaxGroup %in% c("2","4")
        radioMaxGroup <- input$radioMaxGroup
        boxGraph <- input$radioMaxGroup != "4"
      }
      
      return(makePlots(boxData, noLegend, boxGraph, siteToFind))
    })

#############################################################    
    
    output$mymap <- leaflet::renderLeaflet({
      
      leaflet() %>%
        addProviderTiles("Esri.WorldPhysical") %>%
        setView(lng = -83.5, lat = 44.5, zoom=5) 
      
    })
    
    observe({
      
      sumStat <- summaryFile 
      
      if(is.null(input$data)){
        sumStat <- sumStat
      } else if (input$data == "Passive Samples"){
        sumStat <- chemGroup() %>%
          group_by(site) %>%
          summarise(maxEAR = max(EAR)) %>%
          mutate(nSamples = 1) %>%
          mutate(freq = NA)
      }
      
      if(is.null(input$sites) | input$sites == "All"){
        siteToFind <- summaryFile$site
      } else if (input$sites == "Potential 2016"){
        siteToFind <- df2016$Site
        sumStat <- sumStat %>%
          filter(site %in% siteToFind)
        
      } else {
        siteToFind <- input$sites
        sumStat <- sumStat %>%
          filter_(paste0("site == '", siteToFind, "'")) %>%
          data.frame()
      }
      
      mapData <- right_join(stationINFO[,c("shortName", "Station.Name", "dec.lat.va","dec.long.va")], sumStat, by=c("shortName"="site"))
      mapData <- mapData[!is.na(mapData$dec.lat.va),]
      
      col_types <- c("darkblue","dodgerblue","green","yellow","orange","red","brown")
      
      # cols <- colorNumeric(col_types, domain = leg_vals)

      
      if(nrow(mapData) > 1){
        leg_vals <- unique(as.numeric(quantile(mapData$maxEAR, probs=c(0,0.01,0.1,0.25,0.5,0.75,0.9,.99,1), na.rm=TRUE)))
        pal = colorBin(col_types, mapData$maxEAR, bins = leg_vals)
        rad <- 1.5*seq(5000,20000, 1000)
        mapData$sizes <- rad[as.numeric(cut(mapData$nSamples, breaks=16))]
      } else {
        leg_vals <- unique(as.numeric(quantile(c(0,mapData$maxEAR), probs=c(0,0.01,0.1,0.25,0.5,0.75,0.9,.99,1), na.rm=TRUE)))
        pal = colorBin(col_types, c(0,mapData$maxEAR), bins = leg_vals)
        mapData$sizes <- 1.5*12000
      }
      
      
      map <- leafletProxy("mymap", data=mapData) %>%
        clearShapes() %>%
        clearControls() %>%
        addCircles(lat=~dec.lat.va, lng=~dec.long.va, 
                   popup=paste0('<b>',mapData$Station.Name,"</b><br/><table>",
                                "<tr><td>maxEAR: </td><td>",mapData$maxEAR,'</td></tr>',
                                "<tr><td>Number of Samples: </td><td>",mapData$nSamples,'</td></tr>',
                                "<tr><td>Frequency: </td><td>",mapData$freq,'</td></tr></table>') ,
                   fillColor = ~pal(maxEAR), 
                   weight = 1,
                   color = "black",
                   fillOpacity = 0.8, radius = ~sizes, opacity = 0.8) 
      
      if(length(siteToFind) > 1){
        map <- addLegend(map,
          position = 'bottomleft',
          pal=pal,
          values=~maxEAR,
          opacity = 0.8,
          labFormat = labelFormat(transform = function(x) as.integer(x)),
          title = 'Maximum EAR')
        
      }
      
      map
      
    })
    
############################################################# 
    
    output$groupControl <- renderUI({

      if(is.null(input$groupCol)){
        groupCol <- names(endPointInfo)[20]
      } else {
        groupCol <- input$groupCol
      }
      
      siteToFind <- summaryFile$site

      chemGroup <- chemGroup()
      
      statsOfGroup <-  chemGroup %>%
        filter(site %in% siteToFind) %>%
        mutate(category = choices) %>%
        group_by(site, date,category) %>%
        summarise(sumEAR = sum(EAR)) %>%
        group_by(site,category) %>%
        summarise(max = sum(sumEAR > 0.1)) %>%
        data.frame() %>% 
        group_by(category) %>%
        summarise(max=max(max)) %>%
        arrange(desc(max))

     
      ChoicesInGroup <- names(table(endPointInfo[,groupCol]))
      nEndPointsInChoice <- as.character(table(endPointInfo[,groupCol]))
      
      reorderChoices <- setNames(nEndPointsInChoice, ChoicesInGroup)
      
      dropDownHeader <- paste0(statsOfGroup$category," (",reorderChoices[statsOfGroup$category],")")
      
      selectInput("group", label = "Group in annotation (# End Points)",
                  choices = setNames(statsOfGroup$category,dropDownHeader),
                  multiple = FALSE)

    })
    
    output$TableHeader <- renderUI({
      HTML(paste("<h4>", input$group,"-",input$data,": ",input$sites, "</h4>"))
    })
    
    output$mapFooter <- renderUI({
      if(input$data == "Water Sample"){
        HTML("<h5>Size range represents number of collected samples from 1-64</h5>")
      } else {
        HTML("<h5>One sample per site</h5>")
      }
      
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
      
      if(is.null(input$sites) | input$sites == "All"){
        siteToFind <- summaryFile$site
      } else if (input$sites == "Potential 2016"){
        siteToFind <- df2016$Site
      } else {
        siteToFind <- input$sites
      }
      
      if(radio == "2"){
        word <- "chemicals"
      } else if (radio == "3") {
        word <- "classes"
      } else if(radio == "4"){
        word <- "endPoints"
      } else {
        word <- "choices"
      }
      
      if(length(siteToFind) > 1){
        place <- "per site"
      } else {
        place <- ""
      }
      
      if(length(siteToFind) > 1){
        textUI <- paste("<h5>max = Maximum of number of",word,"with hits per site</h5>",
                        "<h5>mean = Mean number of",word,"with hits per site</h5>",
                        "<h5>nSamples = Number of samples per site</h5>")
      } else {
        textUI <- paste("<h5>hits = Number of",word,"with hits </h5>",
                        "<h5>nSamples = Number of samples per site</h5>")
      }
      HTML(textUI)
      
    })

})