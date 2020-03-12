

 
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
      dplyr::filter(SiteID == siteToFind)
  }

  map <-  leaflet::leafletProxy("mymap", data = mapData) %>%
    leaflet::clearMarkers() %>%
    leaflet::addCircleMarkers(lat=~dec_lat, lng=~dec_lon,
                              popup=paste0('<b>',mapData$`Short Name`,"</b><br/><table>",
                                           "<tr><td>",maxEARWords,": </td><td>",sprintf("%.1f",mapData$meanMax),'</td></tr>',
                                           "<tr><td>Number of Samples: </td><td>",mapData$count,'</td></tr>',
                                           '</table>') ,
                              fillColor = ~pal(meanMax),
                              fillOpacity = 0.8,
                              radius = ~sizes,
                              stroke=FALSE,
                              opacity = 0.8) %>%
    leaflet::fitBounds(~min(dec_lon), ~min(dec_lat),
                       ~max(dec_lon), ~max(dec_lat))

  sum_words <- ifelse(sum_logic, "Sum of","Max")

  if(length(siteToFind) > 1){
    map <- map %>% leaflet::clearControls()

    map <- leaflet::addLegend(map,pal = pal,
                              position = 'bottomleft',
                              values=~meanMax,
                              opacity = 0.8,
                              labFormat = leaflet::labelFormat(digits = 2), #transform = function(x) as.integer(x)),
                              title = paste(sum_words,category,"EAR<br>",
                                            maxEARWords,"at site"))
  }

  map
})



output$mapFooter <- renderUI({

  validate(
    need(!is.null(rawData_data$data), "Please select a data set")
  )

  chemical_summary <- chemical_summary()

  nSamples <- dplyr::select(chemical_summary,site,date) %>%
    dplyr::distinct() %>%
    dplyr::group_by(site) %>%
    dplyr::summarize(count = dplyr::n())

  HTML(paste0("<h5>Size range represents number of samples. Ranges from ", min(nSamples$count,na.rm = TRUE),
              " - ",
              max(nSamples$count,na.rm = TRUE),"</h5>"))

})


mapDataINFO <- reactive({
  
  rawData <- rawData()
  
  chem_site <- rawData$chem_site
  chemical_summary <- chemical_summary()
  
  catType = as.numeric(input$radioMaxGroup)
  category <- c("Biological","Chemical","Chemical Class")[catType]
  
  meanEARlogic <- as.logical(input$meanEAR)
  sum_logic <- as.logical(input$sumEAR)
  
  siteToFind <- unique(chemical_summary$site)
  
  mapDataList <- map_tox_data(chemical_summary, 
                              chem_site = chem_site, 
                              category = category,
                              mean_logic = meanEARlogic,
                              sum_logic = sum_logic)
  
  if(length(siteToFind) == 1){
    if(!is.null(latest_map)){
      mapDataList <- latest_map
    } 
  }
  
  latest_map <<- mapDataList
  shinyAce::updateAceEditor(session, editorId = "mapCode_out", value = mapCode() )
  return(mapDataList)
  
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