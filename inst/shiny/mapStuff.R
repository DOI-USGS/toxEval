output$mymap <- leaflet::renderLeaflet({
  
  input$exampleData
  input$data
  
  isolate({
    map <- leaflet() %>%
      addProviderTiles("CartoDB.Positron") %>%
      setView(lng = -83.5, lat = 44.5, zoom=6)    
  })

})

output$mapFooter <- renderUI({
  
  validate(
    need(!is.null(rawData_data$data), "Please select a data set")
  )
  
  chemical_summary <- chemical_summary()
  
  nSamples <- select(chemical_summary,site,date) %>%
    distinct() %>%
    group_by(site) %>%
    summarize(count = n())
  
  HTML(paste0("<h5>Size range represents number of samples. Ranges from ", min(nSamples$count,na.rm = TRUE),
              " - ", 
              max(nSamples$count,na.rm = TRUE),"</h5>"))
  
})


mapDataINFO <- reactive({
  
  rawData <- rawData()
  
  chem_site <- rawData$chem_site
  chemical_summary <- chemical_summary()
  
  siteToFind <- unique(chemical_summary$site)
  
  if(length(siteToFind) == 1){
    mapDataList <- latest_map
  } else {
    catType = as.numeric(input$radioMaxGroup)
    category <- c("Biological","Chemical","Chemical Class")[catType]
    
    meanEARlogic <- as.logical(input$meanEAR)
    sum_logic <- as.logical(input$sumEAR)
    
    mapDataList <- map_tox_data(chemical_summary, 
                              chem_site = chem_site, 
                              category = category,
                              mean_logic = meanEARlogic,
                              sum_logic = sum_logic)
    latest_map <<- mapDataList
  }
  updateAceEditor(session, editorId = "mapCode_out", value = mapCode() )
  return(mapDataList)
  
})

observe({
  
  validate(
    need(!is.null(rawData_data$data), "Please select a data set")
  )

  rawData <- rawData()
  
  chem_site <- rawData$chem_site
  chemical_summary <- chemical_summary()
  
  siteToFind <- unique(chemical_summary$site)
  
  if(length(siteToFind) == 1){
    chem_site <- chem_site[chem_site$SiteID == siteToFind,]
  }
  
  catType = as.numeric(input$radioMaxGroup)
  
  meanEARlogic <- as.logical(input$meanEAR)
  sum_logic <- as.logical(input$sumEAR)
  
  maxEARWords <- ifelse(meanEARlogic,"Mean","Max")
  category <- c("Biological","Chemical","Chemical Class")[catType]
  
  mapDataList <- mapDataINFO()
  
  mapData <- mapDataList$mapData
  pal <- mapDataList$pal
  
  if(length(siteToFind) == 1){

    mapData <- mapData %>%
      filter(SiteID == siteToFind)
  }

  map <- leafletProxy("mymap", data=mapData) %>%
    clearMarkers() %>%
    addCircleMarkers(lat=~dec_lat, lng=~dec_lon,
                     popup=paste0('<b>',mapData$`Short Name`,"</b><br/><table>",
                                  "<tr><td>",maxEARWords,": </td><td>",sprintf("%.1f",mapData$meanMax),'</td></tr>',
                                  "<tr><td>Number of Samples: </td><td>",mapData$count,'</td></tr>',
                                  '</table>') ,
                     fillColor = ~pal(meanMax),
                     fillOpacity = 0.8,
                     radius = ~sizes,
                     stroke=FALSE,
                     opacity = 0.8)
  
    # if(length(siteToFind) > 1){
    #   map <- map %>%
    #     clearControls() %>%
    #     setView(lng = mean(mapData$dec_lon, na.rm = TRUE), 
    #             lat = mean(mapData$dec_lat, na.rm = TRUE), zoom=6) %>%
    #     fitBounds(lng1 = min(mapData$dec_lon, na.rm = TRUE), 
    #               lat1 = min(mapData$dec_lat, na.rm = TRUE), 
    #               lng2 = max(mapData$dec_lon, na.rm = TRUE), 
    #               lat2 = max(mapData$dec_lat, na.rm = TRUE)) 
    # }
  
  sum_words <- ifelse(sum_logic, "Sum of","Max")
  
  if(length(siteToFind) > 1){
    map <- map %>% clearControls() 
    
    map <- addLegend(map,pal = pal,
                     position = 'bottomleft',
                     values=~meanMax,
                     opacity = 0.8,
                     labFormat = labelFormat(digits = 2), #transform = function(x) as.integer(x)),
                     title = paste(sum_words,category,"EAR<br>",
                                   maxEARWords,"at site"))
  }
  
  map
  
})

mapCode <- reactive({
  
  catType = as.numeric(input$radioMaxGroup)
  
  plot_ND = input$plot_ND
  
  category <- c("Biological","Chemical","Chemical Class")[catType]
  
  mapCode <- paste0(rCodeSetup(),"
make_tox_map(chemical_summary, 
        chem_site = tox_list$chem_site,
        category = '",category,"',
        mean_logic = ",as.logical(input$meanEAR),")")

  return(mapCode)
  
})
