library(dplyr)
library(ggplot2)
library(DT)
library(leaflet)
library(data.table)
library(tidyr)
library(toxEval)
library(RColorBrewer)
library(grid)

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

pathToApp <- system.file("extdata", package="toxEval")

summaryFile <- readRDS(file.path(pathToApp,"summary.rds"))
endPoint <- readRDS(file.path(pathToApp,"endPoint.rds"))
chemicalSummaryWS <- readRDS(file.path(pathToApp,"chemicalSummary.rds"))
chemicalSummaryPS <- readRDS(file.path(pathToApp,"chemicalSummaryPassive.rds"))
chemicalSummaryDL <- readRDS(file.path(pathToApp,"chemSummeryDL.rds"))
chemicalSummaryPS$date <- rep(as.POSIXct(as.Date("1970-01-01")),nrow(chemicalSummaryPS))
df2016 <- readRDS(file.path(pathToApp,"df2016.rds"))

choicesPerGroup <- apply(endPointInfo[,-3], 2, function(x) length(unique(x)))
groupChoices <- paste0(names(choicesPerGroup)," (",choicesPerGroup,")")

uniqueClasses <- unique(c(unique(chemicalSummaryPS$class),unique(chemicalSummaryWS$class),unique(chemicalSummaryDL$class)))

#Missing Duluth info:

stationINFO$lakeColor[is.na(stationINFO$lakeColor)] <- "brown"
stationINFO$dec.lat.va[is.na(stationINFO$dec.lat.va)] <- 46.7
stationINFO$dec.long.va[is.na(stationINFO$dec.long.va)] <- -92.4


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

makePlots <- function(boxData, noLegend, boxPlot, uniqueClasses){
  
  siteToFind <- unique(boxData$site)

  if(length(siteToFind) > 1){
    graphData <- boxData %>%
      filter(!is.na(cat)) %>%
      group_by(site,date,cat) %>%
      summarise(sumEAR=sum(EAR)) %>%
      data.frame() %>%
      group_by(site, cat) %>%
      summarise(maxEAR=max(sumEAR),
                meanEAR=mean(sumEAR)) %>%
      data.frame() %>%
      mutate(cat=as.character(cat)) %>%
      gather(stat, value, -site, -cat) %>%
      mutate(cat=factor(cat, levels=uniqueClasses))
    
    orderColsBy <- graphData %>%
      filter(stat == "meanEAR") %>%
      mutate(cat = as.character(cat)) %>%
      group_by(cat) %>%
      summarise(median = quantile(value[value != 0],0.5)) %>%
      filter(!is.na(cat)) %>%
      arrange(median)
    
    orderedLevels <- orderColsBy$cat[!is.na(orderColsBy$median)] #The !is.na takes out any category that was all censo
    orderedLevels <- c(levels(graphData$cat)[!(levels(graphData$cat) %in% orderedLevels)],orderedLevels)
    
    graphData$reorderedCat <- factor(as.character(graphData$cat), levels=orderedLevels)
    
    graphData <- filter(graphData, !is.na(cat))
    
    nlabels <- table(graphData$cat[graphData$stat == "meanEAR"])
    
    lowerPlot <- ggplot(graphData[graphData$stat == "meanEAR",])+
      scale_y_log10("Mean EAR Per Site")
    
    labelsText <<- "nSites"
    
  } else {
    
    siteData <- boxData %>%
      group_by(site,date,cat) %>%
      summarise(sumEAR=sum(EAR)) %>%
      data.frame() %>%
      mutate(cat=factor(cat, levels=uniqueClasses)) %>%
      rename(value=sumEAR)%>%
      filter(!is.na(cat))
    
    orderColsBy <- siteData %>%
      mutate(cat = as.character(cat)) %>%
      group_by(cat) %>%
      summarise(median = median(value[value != 0])) %>%
      filter(!is.na(cat)) %>%
      arrange(median)
    
    orderedLevels <- orderColsBy$cat[!is.na(orderColsBy$median)] #The !is.na takes out any category that was all censo
    orderedLevels <- c(levels(siteData$cat)[!(levels(siteData$cat) %in% orderedLevels)],orderedLevels)
    
    siteData$reorderedCat <- factor(as.character(siteData$cat), levels=orderedLevels)
      
    lowerPlot <- ggplot(siteData)+
      scale_y_log10("Sum of EAR")
    
    nlabels <- table(siteData$cat)
    labelsText <<- "nSamples"
  }
  
  if(!boxPlot){
    lowerPlot <- lowerPlot + geom_point(aes(x=reorderedCat, y=value, color=cat, size=3))
  } else {
    lowerPlot <- lowerPlot + geom_boxplot(aes(x=reorderedCat, y=value, fill=cat)) 
  }
  
  lowerPlot <- lowerPlot + 
    theme(legend.position = "none") +
    theme(axis.text.y = element_text(size=10)) +
    xlab("") +
    geom_hline(yintercept = 0.1, linetype="dashed", color="red")  +
    scale_x_discrete(drop=FALSE) +
    scale_fill_discrete(drop=FALSE)

  ymin <<- 10^(ggplot_build(lowerPlot)$panel$ranges[[1]]$y.range)[1]
  
  namesToPlot <<- names(nlabels)
  nSamples <<- as.character(nlabels)
  
  lowerPlot <- lowerPlot + 
    geom_text(data=data.frame(), aes(x=namesToPlot, y=ymin,label=nSamples),size=3)  +
    geom_text(data=data.frame(), aes(label = c(labelsText,"Hit Threshold"), 
                                     y = c(ymin,0.1), x = c(Inf,Inf),size=3, vjust = -1)) + 
    coord_flip() 
  
  # Code to override clipping
  lowerPlot <- ggplot_gtable(ggplot_build(lowerPlot))
  lowerPlot$layout$clip[lowerPlot$layout$name == "panel"] <- "off"

  siteLimits <- stationINFO %>%
    filter(shortName %in% unique(boxData$site))
  
  if(length(siteToFind) > 1){

    graphData <- graphData %>%
      arrange(as.character(cat))
    
    upperPlot <- ggplot(graphData, aes(x=site, y=value, fill = cat)) +
      geom_bar(stat="identity") +
      facet_wrap(~stat, nrow=2, ncol=1, scales = "free_y") + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.25, 
                                       colour=siteLimits$lakeColor)) +
      scale_x_discrete(limits=siteLimits$Station.shortname,drop=FALSE) +
      xlab("") +
      ylab("EAR") +
      # scale_y_log10("EAR") +
      scale_fill_discrete("", drop=FALSE) 
    
    if(noLegend){
      upperPlot <- upperPlot + 
      guides(fill=FALSE) 
    }

  } else {
    
    chemGroupBPOneSite <- boxData %>%
      select(-site)  %>%
      mutate(cat=factor(cat, levels=uniqueClasses)) %>%
      arrange(as.character(cat)) 
    
    upperPlot <- ggplot(chemGroupBPOneSite, aes(x=date, y=EAR, fill = cat)) +
      geom_bar(stat="identity")+
      theme(axis.text.x=element_blank(),
            axis.ticks=element_blank())+
      xlab("Individual Samples") + 
      ylab("EAR") +
      scale_fill_discrete("", drop=FALSE) +
      labs(fill="") 
      # scale_y_log10("EAR")
    
    if(noLegend){
      upperPlot <- upperPlot +
      guides(fill=FALSE) 
    } 
  }
  
  return(list(upperPlot=upperPlot, lowerPlot=lowerPlot))
}

shinyServer(function(input, output,session) {
  
  choices <- reactive({
    if (input$data == "Duluth" ){
      choices =  c("All",stationINFO$shortName[is.na(stationINFO$queryTime)])
    } else {
      choices =  c("All","Potential 2016",summaryFile$site)
    }
    choices
  })

  observe({
    updateSelectInput(session, "sites", choices = choices())
  })

#############################################################   
    chemicalSummary <- reactive({
      
      if (input$data == "Water Sample"){
        chemicalSummary <- chemicalSummaryWS
      } else if (input$data == "Passive Samples"){
        chemicalSummary <- chemicalSummaryPS
      } else if (input$data == "Duluth"){
        chemicalSummary <- chemicalSummaryDL
      }
      
      chemicalSummary 
      
    })
  
    chemicalSummaryFiltered <- reactive({
      
      chemicalSummary <- chemicalSummary()

      groupCol <- input$groupCol

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
                                   "assay_component_endpoint_name", names(endPointInfo)[4]) %>%
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

      groupCol <- input$groupCol

      radio <- input$radioMaxGroup
      
      statsOfColumn <- chemicalSummary %>%
        rename(assay_component_endpoint_name=endPoint) %>%
        filter(assay_component_endpoint_name %in% endPointInfo$assay_component_endpoint_name ) %>%
        data.table()%>%
        left_join(data.table(endPointInfo[,c("assay_component_endpoint_name", groupCol)]), by = "assay_component_endpoint_name") %>%
        data.frame() %>%
        mutate(site = siteKey[site]) %>%
        rename(endPoint = assay_component_endpoint_name) %>%
        select_("hits","EAR","chnm","class","date","choices"=groupCol,"site","endPoint","endPointValue") 
      
      if (input$sites == "All"){
        siteToFind <- unique(statsOfColumn$site)
      } else if (input$sites == "Potential 2016"){
        siteToFind <- df2016$Site
        statsOfColumn <- statsOfColumn  %>%
          filter(site %in% siteToFind)
      } else {
        siteToFind <- input$sites
        if(any(siteToFind %in% statsOfColumn$site)){
          statsOfColumn <- statsOfColumn  %>%
            filter(site %in% siteToFind)
        } else {
          siteToFind <- unique(statsOfColumn$site)
        }
      }
      
      if(length(siteToFind) == 1){

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
      
      groupCol <- input$groupCol

      chemGroup <- chemicalSummary %>%
        rename(assay_component_endpoint_name=endPoint) %>%
        filter(assay_component_endpoint_name %in% endPointInfo$assay_component_endpoint_name ) %>%
        data.table() %>%
        left_join(data.table(endPointInfo[,c("assay_component_endpoint_name", groupCol)]), by = "assay_component_endpoint_name") %>%
        data.frame() %>%
        mutate(site = siteKey[site]) %>%
        rename(endPoint=assay_component_endpoint_name) %>%
        select_("hits","EAR","chnm","class","site","date","choices"=groupCol,"endPoint", "endPointValue")
      
    })
    
    statsOfGroupOrdered <- reactive({
      
      groupCol <- input$groupCol

      radioMaxGroup <- input$radioMaxGroup
      
      chemGroup <- chemGroup()   
      
      if (input$sites == "All"){
        siteToFind <- unique(chemGroup$site)
        statsOfGroup <-  chemGroup
      } else if (input$sites == "Potential 2016"){
        siteToFind <- df2016$Site
        statsOfGroup <-  chemGroup %>%
          filter(site %in% siteToFind)
      } else {
        if(any(input$sites %in% chemGroup$site)){
          siteToFind <- input$sites
          statsOfGroup <-  chemGroup %>%
            filter(site %in% siteToFind)
        } else {
          siteToFind <- unique(chemGroup$site)
          statsOfGroup <-  chemGroup
        }
      }
      
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

      radio <- input$radio

      groupCol <- input$groupCol
      
      group <- input$group
      
      chemicalSummaryFiltered <- chemicalSummaryFiltered()
      
      if (input$sites == "All"){
        siteToFind <- unique(chemicalSummaryFiltered$site)
      } else if (input$sites == "Potential 2016"){
        siteToFind <- df2016$Site
        
        chemicalSummaryFiltered <- chemicalSummaryFiltered %>%
          filter(site %in% siteToFind) 
      } else {
        if(any(input$sites %in% chemicalSummaryFiltered$site)){
          siteToFind <- input$sites
          
          chemicalSummaryFiltered <- chemicalSummaryFiltered %>%
            filter(site %in% siteToFind) 
        } else {
          siteToFind <- unique(chemicalSummaryFiltered$site)
        }
        
      }
      
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
                                      options = list(dom = 'ft',
                                                     pageLength = nrow(statTable), 
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
  
      if (input$sites == "All"){
        siteToFind <- unique(statsOfGroupOrdered$site)
      } else if (input$sites == "Potential 2016"){
        siteToFind <- df2016$Site
      } else {
        if(any(input$sites %in% statsOfGroupOrdered$site)){
          siteToFind <- input$sites
          statsOfGroupOrdered <- statsOfGroupOrdered %>%
            filter(siteToFind %in% site)
        } else {
          siteToFind <- unique(statsOfGroupOrdered$site)
        }
      }

      if(length(siteToFind) > 1){
        
        colToSort <- 1
        if (input$data == "Water Sample"){
          
          summaryFile <- filter(summaryFile, site %in% siteToFind)
          statsOfGroupOrdered <- left_join(summaryFile[,c("site","nSamples")], statsOfGroupOrdered, by="site")
          
          colToSort <- 2
        }
        
        meanChem <- grep("mean",names(statsOfGroupOrdered))
        maxChem <- grep("max",names(statsOfGroupOrdered))

        tableGroup <- DT::datatable(statsOfGroupOrdered, 
                                      rownames = FALSE,
                                      filter = 'top',
                                      options = list(dom = 'ft',
                                                     pageLength = nrow(statsOfGroupOrdered), 
                                                     order=list(list(colToSort,'desc'))))

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
      
      colToSort <- 1
      if("nSamples" %in% names(statCol)){
        colToSort <- 2
      }
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
                                 options = list(dom = 'ft',
                                                pageLength = nrow(statCol),
                                                order=list(list(colToSort,'desc'))))

      
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
      print(grid.draw(plots()$lowerPlot))
    })
    
    output$graphGroup <- renderPlot({ 
      print(grid.draw(groupPlots()$lowerPlot))
    })
    
    boxData2 <- reactive({
      
      radio <- input$radio

      groupCol <- input$groupCol
      
      group <- input$group
      
      chemGroup <- chemGroup()
      
      if (input$sites == "All"){
        siteToFind <- unique(chemGroup$site)
      } else if (input$sites == "Potential 2016"){
        siteToFind <- df2016$Site
        chemGroup <- chemGroup %>%
          filter(site %in% siteToFind) 
      } else {
        if(any(input$sites %in% chemGroup$site)){
          siteToFind <- input$sites
          chemGroup <- chemGroup %>%
            filter(site %in% siteToFind) 
        } else {
          siteToFind <- unique(chemGroup$site)
        }
      }
      
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
      
      boxData2 <- chemGroupBP_group %>%
        filter(choices %in% group) 
        # filter(EAR > 0)
      
      return(boxData2)
      
    })
    
    plots <- reactive({
      
      boxData <- boxData2()
      
      siteToFind <- unique(boxData$site)
      
      noLegend <- TRUE

      if(input$radio == "2"){
        levelsToGraph <- uniqueClasses 
      } else {
        levelsToGraph <- unique(boxData$cat)
      }
      
      boxPlot <- !(input$data == "Passive Samples" & length(siteToFind) == 1)
      
      return(makePlots(boxData, noLegend, boxPlot, levelsToGraph))
      
    })

    boxData <- reactive({
      
      radioMaxGroup <- input$radioMaxGroup
      
      chemGroup <- chemGroup()

      if (input$sites == "All"){
        siteToFind <- unique(chemGroup$site)
      } else if (input$sites == "Potential 2016"){
        siteToFind <- df2016$Site
        chemGroup <- chemGroup %>%
          filter(site %in% siteToFind)
      } else {
        if(any(input$sites %in% chemGroup$site)){
          siteToFind <- input$sites
          chemGroup <- chemGroup %>%
            filter(site %in% siteToFind)
        } else {
          siteToFind <- unique(chemGroup$site)
        }
      }
            
      boxData <- chemGroup %>%
        mutate(date = factor(as.numeric(date) - min(as.numeric(date))))
      
      if(radioMaxGroup == "1"){
        boxData <- mutate(boxData, cat=choices)
      } else if (radioMaxGroup == "2"){
        boxData <- mutate(boxData, cat=chnm)
      } else {
        boxData <- mutate(boxData, cat=class)
      }

      
      return(boxData)
      
    })
    
    groupPlots <- reactive({
      
      boxData <- boxData()
      
      noLegend <- TRUE
      boxGraph <- TRUE
      levelsToGraph <- unique(boxData$cat)[!is.na(unique(boxData$cat))]
      siteToFind <- unique(boxData$site)

      radioMaxGroup <- input$radioMaxGroup
      boxGraph <-  !(input$data == "Passive Samples" & length(siteToFind) == 1) 
      if(input$radioMaxGroup == "3"){
        levelsToGraph <- uniqueClasses
      }
      
      return(makePlots(boxData, noLegend, boxGraph, levelsToGraph))
    })

#############################################################    
    
    output$mymap <- leaflet::renderLeaflet({
      
      leaflet(height = "50px") %>%
        # addProviderTiles("Esri.WorldPhysical") %>%
        addProviderTiles("CartoDB.Positron") %>%
        setView(lng = -83.5, lat = 44.5, zoom=6) 
      
    })
    
    observe({
      
      sumStat <- summaryFile 
      
      if (input$data == "Passive Samples"){
        sumStat <- chemGroup() %>%
          group_by(site) %>%
          summarise(maxEAR = max(EAR)) %>%
          mutate(nSamples = 1) %>%
          mutate(freq = NA)        
      } else if (input$data == "Duluth"){
        sumStat <- chemGroup() %>%
          group_by(site, date) %>%
          summarise(maxEAR = max(EAR),
                    hits=as.numeric(any(hits > 0))) %>%
          data.frame() %>%
          group_by(site) %>%
          summarise(nSamples = n(),
                    maxEAR=max(maxEAR,na.rm=TRUE),
                    freq=sum(hits)/n()) %>%
          data.frame()
      } else {
        sumStat <- sumStat
      }
      
      if (input$sites == "All"){
        siteToFind <- unique(sumStat$site)
      } else if (input$sites == "Potential 2016"){
        siteToFind <- df2016$Site
        sumStat <- sumStat %>%
          filter(site %in% siteToFind)
        
      } else {
        if(any(input$sites %in%sumStat$site)){
          siteToFind <- input$sites
          sumStat <- sumStat %>%
            filter(site %in% siteToFind)
        } else {
          siteToFind <- unique(sumStat$site)
        }

      }
      
      mapData <- right_join(stationINFO[,c("shortName", "Station.Name", "dec.lat.va","dec.long.va")], sumStat, by=c("shortName"="site"))
      mapData <- mapData[!is.na(mapData$dec.lat.va),]
      
      col_types <- c("darkblue","dodgerblue","chartreuse2","gold1","orange","brown","red")

      if(nrow(mapData) > 1){
        leg_vals <- unique(as.numeric(quantile(mapData$maxEAR, probs=c(0,0.01,0.1,0.25,0.5,0.75,0.9,.99,1), na.rm=TRUE)))
        pal = colorBin(col_types, mapData$maxEAR, bins = leg_vals)
        rad <-3*seq(1,4,length.out = 16)
        # rad <- 1.5*seq(5000,20000, 1000)
        mapData$sizes <- rad[as.numeric(cut(mapData$nSamples, breaks=16))]
      } else {
        leg_vals <- unique(as.numeric(quantile(c(0,mapData$maxEAR), probs=c(0,0.01,0.1,0.25,0.5,0.75,0.9,.99,1), na.rm=TRUE)))
        pal = colorBin(col_types, c(0,mapData$maxEAR), bins = leg_vals)
        mapData$sizes <- 3
        # mapData$sizes <- 1.5*12000
      }
      
      
      map <- leafletProxy("mymap", data=mapData) %>%
        # clearShapes() %>%
        clearMarkers() %>%
        clearControls() %>%
        addCircleMarkers(lat=~dec.lat.va, lng=~dec.long.va, 
                   popup=paste0('<b>',mapData$Station.Name,"</b><br/><table>",
                                "<tr><td>maxEAR: </td><td>",sprintf("%.1f",mapData$maxEAR),'</td></tr>',
                                "<tr><td>Number of Samples: </td><td>",mapData$nSamples,'</td></tr>',
                                "<tr><td>Frequency: </td><td>",sprintf("%.1f",mapData$freq),'</td></tr></table>') ,
                   fillColor = ~pal(maxEAR), 
                   # weight = 1,
                   # color = "black",
                   fillOpacity = 0.8, 
                   radius = ~sizes, 
                   stroke=FALSE,
                   opacity = 0.8) 
      
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
    
    observe({

     groupCol <- input$groupCol

     chemGroup <- chemGroup()
      
     orderBy <- chemGroup %>%
        group_by(choices) %>%
        summarise(max = max(EAR),
                  nEndPoint = length(unique(endPoint))) %>%
        arrange(desc(max)) %>%
        filter(!is.na(choices))
      
     dropDownHeader <- paste0(orderBy$choices," (",orderBy$nEndPoint,")")
      
     updateSelectInput(session, "group",
      choices = setNames(orderBy$choices,dropDownHeader)
     )

    })
    
    output$TableHeader <- renderUI({
      HTML(paste("<h4>", input$group,"-",input$data,": ",input$sites, "</h4>"))
    })
    
    output$mapFooter <- renderUI({
      if(input$data == "Water Sample"){
        HTML("<h5>Size range represents number of collected samples from 1-64</h5>")
      } else if (input$data == "Duluth") {
        HTML("<h5>Size range represents number of collected samples from 1-26</h5>")
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

      radio <- input$radioMaxGroup
      
      if(input$sites == "All"){
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
    
    output$endpointGraph <- renderPlot({ 

      boxData <- boxData()
    
      filterCat <- input$filterCat
      
      boxData3 <- boxData %>%
        filter(cat == filterCat) %>%
        filter(endPoint %in% unique(endPoint[EAR > 0.1]))

      stackedPlot <- ggplot(boxData3, aes(x = cat, y = EAR))+
        scale_y_log10("Mean EAR Per Site") +
        geom_boxplot(aes(x=endPoint, y=EAR)) +
        coord_flip() +
        xlab("") +
        theme(axis.text.y = element_text(vjust = .25,hjust=1)) +
        geom_hline(yintercept = 0.1, linetype="dashed", color="red")
      
      print(stackedPlot)

      
    })
    
    output$endpointGraph2 <- renderPlot({ 
      
      boxData2 <- boxData2()

      filterCat2 <- input$filterCat2
      
      boxData4 <- boxData2 %>%
        filter(cat == filterCat2) %>%
        filter(endPoint %in% unique(endPoint[EAR > 0.1]))
      
      stackedPlot <- ggplot(boxData4)+
        scale_y_log10("Mean EAR Per Site") +
        geom_boxplot(aes(x=endPoint, y=EAR)) + 
        coord_flip() +
        xlab("") +
        theme(axis.text.y = element_text(vjust = .25,hjust=1))+
        geom_hline(yintercept = 0.1, linetype="dashed", color="red")
      
      print(stackedPlot)
      
    })
    
    output$endpointTable2 <- DT::renderDataTable({    
      boxData2 <- boxData2()

      filterCat2 <- input$filterCat2
      
      boxData4 <- boxData2 %>%
        filter(cat == filterCat2) %>%
        # filter(endPoint %in% unique(endPoint[EAR > 0.1])) %>%
      group_by(site, date, endPoint) %>%
        summarise(sumEAR = sum(EAR),
                  nHits = sum(hits)) %>%
        group_by(site, endPoint) %>%
        summarise(maxEAR = max(sumEAR),
                  meanEAR = mean(sumEAR),
                  sumHits = sum(nHits),
                  freq = sum(nHits > 0)/n()) %>%
        data.frame()
      
      if(length(unique(boxData4$site)) > 1){
        boxData4 <- group_by(boxData4, endPoint) %>%
          summarise(nHits = sum(sumHits > 0)) 
      } else {
        boxData4 <- boxData4 %>%
          select(endPoint,freq,maxEAR)
      }
      
      tableGroup <- DT::datatable(boxData4, 
                                  rownames = FALSE,
                                  filter = 'top',
                                  options = list(pageLength = nrow(boxData4), 
                                                 order=list(list(1,'desc'))))
      
      tableGroup <- formatStyle(tableGroup, names(boxData4)[2], 
                                background = styleColorBar(range(boxData4[,2],na.rm=TRUE), 'goldenrod'),
                                backgroundSize = '100% 90%',
                                backgroundRepeat = 'no-repeat',
                                backgroundPosition = 'center' ) 
      tableGroup
      
    })
    
    output$hitsTable <- DT::renderDataTable({    
      
      boxData <- boxData()
      
      tableData <- boxData %>%
        group_by(site, choices, class) %>%
        summarize(hits = any(hits > 0)) %>%
        group_by(choices, class) %>%
        summarize(nSites = sum(hits)) %>%
        data.frame() %>%
        reshape(idvar="choices",timevar="class", direction="wide") 
      
      names(tableData) <- gsub("nSites.","",names(tableData))
      names(tableData)[1] <- "Groups"
      
      sumOfColumns <- colSums(tableData[-1],na.rm = TRUE)
      orderData <- order(sumOfColumns,decreasing = TRUE) 
      orderData <- orderData[sumOfColumns[orderData] != 0] + 1
      
      tableData <- tableData[,c(1,orderData)]
      colors <- brewer.pal(9,"Blues") #"RdYlBu"
      
      groups <- tableData$Groups
      
      tableData <- tableData[!is.na(groups),-1]
      rownames(tableData) <- groups[!is.na(groups)]

      cuts <- seq(0,max(as.matrix(tableData),na.rm=TRUE),length.out = 8)
            
      tableData1 <- DT::datatable(tableData, 
                                  rownames = TRUE,
                                  options = list(dom = 't',
                                                 pageLength = nrow(tableData), 
                                                 order=list(list(1,'desc'))))

      for(i in 1:ncol(tableData)){
        tableData1 <- formatStyle(tableData1, columns = names(tableData)[i], 
                    backgroundColor = styleInterval(cuts = cuts,values = colors),
                    color = styleInterval(40,values = c("black","white")),
                    `font-size` = '15px')
      }


      tableData1
      
    })

    observe({
      boxData <- boxData()%>%
        filter(EAR >= 0.1)
      
      updateSelectInput(session, "filterCat",
            choices = unique(boxData$cat)) 
    })
    
    observe({
      boxData2 <- boxData2() %>%
        filter(EAR >= 0.1)
      
      updateSelectInput(session, "filterCat2",
             choices = unique(boxData2$cat)
      )
      
    })

})