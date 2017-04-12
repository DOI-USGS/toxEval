output$mymap <- leaflet::renderLeaflet({
  
  map <- leaflet(height = "50px") %>%
    addProviderTiles("CartoDB.Positron") %>%
    setView(lng = -83.5, lat = 44.5, zoom=6)
  
  map
  
})

output$mapFooter <- renderUI({
  
  validate(
    need(!is.null(input$data), "Please select a data set")
  )
  
  chemicalSummary <- chemicalSummary()
  catType = as.numeric(input$radioMaxGroup)
  meanEARlogic <- as.logical(input$meanEAR)  
  hit_threshold <- hitThresValue()
  
  statsOfGroupOrdered <- statsOfGroup(chemicalSummary = chemicalSummary,
                                      category = c("Biological","Chemical","Chemical Class")[catType],
                                      hit_threshold = hit_threshold)
  statsOfGroupOrdered <- statsOfGroupOrdered %>%
    group_by(site) %>%
    summarize(nSamples = median(nSamples, na.rm = TRUE))
  
  if(input$radioMaxGroup == "1"){
    word <- "groups"
  } else if (input$radioMaxGroup == "2"){
    word <- "chemicals"
  } else {
    word <- "classes"
  }
  
  HTML(paste0("<h5>Size range represents number of ",word,
              " with hits. Ranges from ", min(statsOfGroupOrdered$nSamples,na.rm = TRUE),
              " - ", 
              max(statsOfGroupOrdered$nSamples,na.rm = TRUE),"</h5>"))
  
})


observe({
  
  validate(
    need(!is.null(input$data), "Please select a data set")
  )
  
  rawData <- rawData()
  chem_site <- rawData$chem_site
  
  chemicalSummary <- chemicalSummary()
  catType = as.numeric(input$radioMaxGroup)
  siteToFind <- unique(chemicalSummary$site)
  
  meanEARlogic <- as.logical(input$meanEAR)
  
  maxEARWords <- ifelse(meanEARlogic,"meanEAR","maxEAR")
  
  mapDataList <- getMapInfo(chemicalSummary, 
                     chem_site = chem_site, 
                     category = c("Biological","Chemical","Chemical Class")[catType],
                     mean_logic = meanEARlogic,
                     hit_threshold = hitThresValue())
  
  mapData <- mapDataList$mapData
  pal <- mapDataList$pal
  
  if(length(siteToFind) == 1){

    mapData <- filter(chem_site, SiteID == siteToFind) %>%
      mutate(nSamples = median(mapData$nSamples),
             meanMax = median(mapData$meanMax),
             sizes = median(mapData$sizes))
  }

  map <- leafletProxy("mymap", data=mapData) %>%
    clearMarkers() %>%
    clearControls() %>%
    addCircleMarkers(lat=~dec_lat, lng=~dec_lon,
                     popup=paste0('<b>',mapData$`Short Name`,"</b><br/><table>",
                                  "<tr><td>",maxEARWords,": </td><td>",sprintf("%.1f",mapData$meanMax),'</td></tr>',
                                  "<tr><td>Number of Samples: </td><td>",mapData$nSamples,'</td></tr>',
                                  # "<tr><td>Frequency: </td><td>",sprintf("%.1f",mapData$freq),'</td></tr>',
                                  # "<tr><td>Number of ",typeWords," with hits: </td><td>",counts,'</td></tr>',
                                  '</table>') ,
                     fillColor = ~pal(meanMax),
                     # weight = 1,
                     # color = "black",
                     fillOpacity = 0.8,
                     radius = ~sizes,
                     stroke=FALSE,
                     opacity = 0.8)
  
  
  if(length(siteToFind) > 1){
    map <- addLegend(map,pal = pal,
                     position = 'bottomleft',
                     values=~meanMax,
                     opacity = 0.8,
                     labFormat = labelFormat(digits = 2), #transform = function(x) as.integer(x)),
                     title = ifelse(meanEARlogic,'Mean EAR','Max EAR'))
  }
  
  map
  
})
