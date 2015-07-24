library(dplyr)
library(magrittr)
library(ggplot2)
library(DT)
library(leaflet)
library(toxEval)
library(data.table)

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
  mutate(lakeColor = c("red","black","green","brown","brown","blue","blue")[as.numeric(lakeCat)] )

endPointInfo <- endPointInfo

pathToApp <- system.file("extdata", package="toxEval")

summary <- readRDS(file.path(pathToApp,"summary.rds"))
endPoint <- readRDS(file.path(pathToApp,"endPoint.rds"))
chemicalSummary <- readRDS(file.path(pathToApp,"chemicalSummary.rds"))

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
      select_("hits","EAR","endPoint",
              chemicalSummary.site,chemicalSummary.Date,groupCol) %>%
      group_by_(chemicalSummary.site, chemicalSummary.Date) %>%
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
    mutate(site = siteKey[site]) %>%
    arrange(desc(maxEAR))
    
  })
  
  output$groupControl <- renderUI({
    
    selectInput("group", label = "Group in column",
                choices = unique(endPointInfo[,input$groupCol])[!is.na(unique(endPointInfo[,input$groupCol]))],
                selected = unique(endPointInfo[,20])[3],
                multiple = FALSE)
  })
  
  output$groupSumTable <- DT::renderDataTable({
#     chemicalSummaryFiltered() %>%
#       group_by(site) %>%
#       summarise(nChem = length(unique(chnm)),
#                 nEndPoints = length(unique(endPoint))) %>%
#       select(nChem,nEndPoints) %>%
#       mutate(nChem = median(nChem),
#              nEndPoints = median(nEndPoints)) %>%
#       distinct()
    data.frame(nChem=39, nEndPoints=538)
      
  }, options = list(searching = FALSE, paging = FALSE))
  
  output$table <- DT::renderDataTable({
    statsOfSum()
  }, options = list(pageLength = 10))
  
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
      mutate(site = factor(site, levels=siteLimits$Station.shortname))
    
    sToxWS <- ggplot(endPointSummBP, aes(x=site, y=sumEAR)) +
      geom_boxplot() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.25, colour=siteLimits$lakeColor), 
            legend.position = "none")+
      scale_x_discrete(limits=siteLimits$Station.shortname) +
      # scale_y_log10(limits=c(0.03,5000)) +
      geom_text(data=data.frame(), aes(x=c(5, 18,31,42,54),
                                       y=-.5, label=c("Superior","Michigan","Huron","Erie","Ontario")),
                colour=factor(c("red","black","green","brown","blue"),
                              levels=c("red","black","green","brown","blue")), size=3)
    print(sToxWS)
  })
  
  output$mymap <- leaflet::renderLeaflet({
    
    leaflet() %>%
      addProviderTiles("Esri.WorldPhysical") %>%
      setView(lng = -83.5, lat = 44.5, zoom=6) 
    
  })
  
  observe({
    
    mapData <- left_join(stationINFO, summary, by=c("shortName"="site"))
    mapData <- mapData[!is.na(mapData$dec.lat.va),]
    mapData <- mapData[!is.na(mapData$nChem),]
    
    if(input$sites != "All"){
      mapData <- mapData[mapData$shortName == input$sites,]
    }
    
    col_types <- c("darkblue","dodgerblue","green","yellow","orange","red")
    leg_vals <- c(0,1,5,10,100,1000,2500)#seq(0,2500, 500) # This is maxEAR 
    
    cols <- colorNumeric(col_types, domain = leg_vals)
    rad <- seq(5000,15000, 500)
    mapData$sizes <- rad[as.numeric(cut(mapData$nChem, c(-1:19)))]
    
    pal = colorBin(col_types, mapData$maxEAR, bins = leg_vals)
    
    leafletProxy("mymap", data=mapData) %>%
      clearShapes() %>%
      clearControls() %>%
      addCircles(lat=~dec.lat.va, lng=~dec.long.va, 
                 popup=
                   paste0('<b>',mapData$Station.Name,"</b><br/><table>",
                          "<tr><td>Frequency</td><td>",sprintf("%1.1f",mapData$freq),'</td></tr>',
                          "<tr><td>nSample</td><td>",mapData$nSamples,'</td></tr>',
                          "<tr><td>maxEAR</td><td>",sprintf("%1.1f",mapData$maxEAR),'</td></tr>',
                          "<tr><td>nChem</td><td>",mapData$nChem,'</td></tr>',
                          "<tr><td>nEndPoints</td><td>",mapData$nEndPoints,'</td></tr></table>')
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