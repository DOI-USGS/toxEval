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

endPointInfo <- endPointInfo
wData <- wData
pCodeInfo <- pCodeInfo

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

endPointInfo <- endPointInfo

pathToApp <- system.file("extdata", package="toxEval")

summary <- readRDS(file.path(pathToApp,"summary.rds"))
endPoint <- readRDS(file.path(pathToApp,"endPoint.rds"))
chemicalSummary <- readRDS(file.path(pathToApp,"chemicalSummary.rds"))

choicesPerGroup <- apply(endPointInfo[,-3], 2, function(x) length(unique(x)))
groupChoices <- paste0(names(choicesPerGroup)," (",choicesPerGroup,")")

shinyServer(function(input, output) {
  
  # chemicalSummaryFiltered <- eventReactive(input$calculate, {
  chemicalSummaryFiltered <- reactive({
    
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
    
  })
  
  # endpointSummary <- eventReactive(input$calculate, {
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
      select_("hits","EAR","endPoint","site","date",groupCol) %>%
      group_by(site, date) %>%
      summarise(sumEAR = sum(EAR),
                nHits = sum(hits)) 
  })
  
  # statsOfSum <- eventReactive(input$calculate, {
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
  
  statsOfColumn <- reactive({
    if(is.null(input$groupCol)){
      groupCol <- names(endPointInfo)[20]
    } else {
      groupCol <- input$groupCol
    }
    
    statsOfColumn <- chemicalSummary %>%
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
  })
  
  output$groupControl <- renderUI({
    
    ChoicesInGroup <- names(table(endPointInfo[,input$groupCol]))

    nEndPointsInChoice <- as.character(table(endPointInfo[,input$groupCol]))
    
    dropDownHeader <- paste0(ChoicesInGroup," (",nEndPointsInChoice,")")
    
    selectInput("group", label = "Group in column (# End Points)",
                choices = setNames(ChoicesInGroup,dropDownHeader),
                selected = unique(endPointInfo[,20])[3],
                multiple = FALSE)
  })

  output$table <- DT::renderDataTable({
    statsOfSumDF <- DT::datatable(statsOfSum(), 
                                rownames = FALSE,
                                colnames = c('Maximum EAR' = 3, 'Sum of Hits' = 4),
                                filter = 'top',
                                options = list(pageLength = 10, 
                                               order=list(list(2,'desc')))) %>%
        formatRound(c("Maximum EAR","meanEAR"), 1) %>%
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
  })
  
  output$tableSumm <- DT::renderDataTable({
    statCol <- statsOfColumn()
    
    maxEARS <- grep("maxEAR",names(statCol))
    
    orderList <- list()
    for(i in 1:length(maxEARS)){
      orderList <- append(orderList, list(maxEARS[i],'desc'))
    }
    
    
    DT::datatable(statCol, 
                  rownames = FALSE,
                  # filter = 'top',
                  options = list(pageLength = 10, 
                                 order=list(orderList)))%>%
      formatRound(names(statCol)[maxEARS], 1) 
  })
  
  output$graph <- renderPlot({
    
    endpointSummary <- endpointSummary()
    
    if(nrow(endpointSummary) == 0){
      endpointSummary$site <- stationINFO$fullSiteID
    }
      
    siteLimits <- stationINFO %>%
      filter(fullSiteID %in% unique(endpointSummary$site))

    
    endPointSummBP <- endpointSummary %>%
      data.frame()%>%
      mutate(site = siteKey[site]) %>%
      mutate(site = factor(site, levels=siteLimits$Station.shortname)) %>%
      mutate(sumEARnoZero = sumEAR) 
    
    ndLevel <- 0.1*min(endPointSummBP$sumEARnoZero[endPointSummBP$sumEARnoZero != 0])
    
    endPointSummBP$sumEARnoZero[endPointSummBP$sumEARnoZero == 0] <- ndLevel
    
    sToxWS <- ggplot(endPointSummBP, aes(x=site, y=sumEARnoZero)) +
      geom_boxplot() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.25, colour=siteLimits$lakeColor), 
            legend.position = "none")+
      scale_x_discrete(limits=siteLimits$Station.shortname) +
      scale_y_log10("sumEAR") +
      coord_cartesian(ylim = c(1, 1.1*max(endPointSummBP$sumEARnoZero))) +
      # ylab("sumEAR") + 
      xlab("") +
      annotation_custom(xmin=-1,xmax=-1,
                        ymin=log10(ndLevel), ymax=log10(ndLevel),
                        grob=textGrob("ND", gp=gpar(fontsize=9), vjust = 0.25)) +
      geom_text(data=data.frame(), 
                aes(x=c(5, 18,31,42,54),y=rep(1.1,5),
                    label=c("Superior","Michigan","Huron","Erie","Ontario")),                
                colour=factor(c("tomato3","black","springgreen3","brown","blue"),
                              levels=c("tomato3","black","springgreen3","brown","blue")), size=3)            
                              
    print(sToxWS)
#     for(i in 1:5){
#       sToxWS <- sToxWS +       
#         annotation_custom(xmin=c(5, 18,31,42,54)[i],xmax=c(5, 18,31,42,54)[i],
#                          ymin=log10(.98), ymax=log10(.98),
#                          grob=textGrob(c("Superior","Michigan","Huron","Erie","Ontario")[i],
#                                        gp=gpar(col=c("red","black","green","brown","blue")[i], fontsize=9)))
# 
#     }

#     g_sToxWS <- ggplotGrob(sToxWS)
#     g_sToxWS$layout$clip[g_sToxWS$layout$name=="panel"] <- "off"
#     grid.draw(g_sToxWS)
    
  })
  
  output$mymap <- leaflet::renderLeaflet({
    
    leaflet() %>%
      addProviderTiles("Esri.WorldPhysical") %>%
      setView(lng = -83.5, lat = 44.5, zoom=6) 
    
  })
  
  observe({

    sumStat <- statsOfSum()
    
    if(nrow(sumStat) == 0){
      sumStat <- summary
      sumStat$sumHits <- sumStat$nChem
    }
    
    mapData <- right_join(stationINFO, sumStat, by=c("shortName"="site"))
    mapData <- mapData[!is.na(mapData$dec.lat.va),]
    # mapData <- mapData[!is.na(mapData$nChem),]
    
    if(input$sites != "All"){
      mapData <- mapData[mapData$shortName == input$sites,]
    }
    
    col_types <- c("darkblue","dodgerblue","green","yellow","orange","red","brown")
    leg_vals <- unique(as.numeric(quantile(mapData$maxEAR, probs=c(0,0.01,0.1,0.25,0.5,0.75,0.9,.99,1), na.rm=TRUE)))
    #c(0,1,5,10,100,1000,2500)#seq(0,2500, 500) # This is maxEAR 
    
    cols <- colorNumeric(col_types, domain = leg_vals)
    rad <- seq(5000,15000, 500)
    # mapData$sizes <- rad[as.numeric(cut(mapData$nChem, c(-1:19)))]
    mapData$sizes <- rad[as.numeric(
                              cut(mapData$sumHits, 
                                  c(-1,unique(as.numeric(
                                    quantile(mapData$sumHits, probs=c(0.01,0.1,0.25,0.5,0.75,0.9,.99), na.rm=TRUE)
                                    )),(max(mapData$sumHits,na.rm=TRUE)+1))
                                  ))]
    pal = colorBin(col_types, mapData$maxEAR, bins = leg_vals)
    
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