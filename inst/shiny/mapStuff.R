output$mymap <- leaflet::renderLeaflet({
  
  map <- leaflet(height = "1000px") %>%
    addProviderTiles("CartoDB.Positron") %>%
    setView(lng = -83.5, lat = 44.5, zoom=6)
  
  map
  
})

output$mapFooter <- renderUI({
  
  validate(
    need(!is.null(input$data), "Please select a data set")
  )
  
  chemicalSummary <- chemicalSummary()
  
  nSamples <- select(chemicalSummary,site,date) %>%
    distinct() %>%
    group_by(site) %>%
    summarize(count = n())
  
  HTML(paste0("<h5>Size range represents number of samples. Ranges from ", min(nSamples$count,na.rm = TRUE),
              " - ", 
              max(nSamples$count,na.rm = TRUE),"</h5>"))
  
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
                     mean_logic = meanEARlogic)
  
  mapData <- mapDataList$mapData
  pal <- mapDataList$pal
  
  if(length(siteToFind) == 1){

    mapData <- filter(chem_site, SiteID == siteToFind) %>%
      mutate(nSamples = median(mapData$count),
             meanMax = median(mapData$meanMax),
             sizes = median(mapData$sizes))
  }

  map <- leafletProxy("mymap", data=mapData) %>%
    clearMarkers() %>%
    clearControls() %>%
    setView(lng = mean(mapData$dec_lon, na.rm = TRUE), 
            lat = mean(mapData$dec_lat, na.rm = TRUE), zoom=6) %>%
    fitBounds(lng1 = min(mapData$dec_lon, na.rm = TRUE), 
              lat1 = min(mapData$dec_lat, na.rm = TRUE), 
              lng2 = max(mapData$dec_lon, na.rm = TRUE), 
              lat2 = max(mapData$dec_lat, na.rm = TRUE)) %>%
    addCircleMarkers(lat=~dec_lat, lng=~dec_lon,
                     popup=paste0('<b>',mapData$`Short Name`,"</b><br/><table>",
                                  "<tr><td>",maxEARWords,": </td><td>",sprintf("%.1f",mapData$meanMax),'</td></tr>',
                                  "<tr><td>Number of Samples: </td><td>",mapData$count,'</td></tr>',
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
