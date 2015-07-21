library(dplyr)
library(magrittr)
library(ggplot2)
library(DT)
library(leaflet)
library(toxEval)

endPointInfo <- endPointInfo
wData <- wData
pCodeInfo <- pCodeInfo

packagePath <- system.file("extdata", package="toxEval")
filePath <- file.path(packagePath, "stationINFO.RData")
load(file=filePath)
siteKey <- setNames(stationINFO$shortName, stationINFO$fullSiteID)

endPointInfo <- endPointInfo

pathToApp <- system.file("extdata", package="toxEval")

summary <- readRDS(file.path(pathToApp,"summary.rds"))
endPoint <- readRDS(file.path(pathToApp,"endPoint.rds"))
chemicalSummary <- readRDS(file.path(pathToApp,"chemicalSummary.rds"))

shinyServer(function(input, output) {
  
  endpointSummary <- reactive({
    
    if(is.null(input$groupCol)){
      groupCol <- names(endPointInfo)[10]
    } else {
      groupCol <- input$groupCol
    }
    
    if(is.null(input$group)){
      group <- unique(endPointInfo[,10])[5]
    } else {
      group <- input$group
    }
    
    endPointInfoSub <- select_(endPointInfo, "assay_component_endpoint_name", groupCol) %>%
      filter_(paste0(groupCol," == '", group, "'")) %>%
      unique()
    
    endpointSummary <- chemicalSummary %>%
      filter(endPoint %in% endPointInfoSub$assay_component_endpoint_name ) %>%
      left_join(endPointInfoSub, by=c("endPoint"="assay_component_endpoint_name")) %>%
      select_("hits","EAR","endPoint","site","date",groupCol) %>%
      group_by(site,date) %>%
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
    mutate(site = siteKey[site]) %>%
    arrange(desc(maxEAR))
    
  })
  
  output$groupControl <- renderUI({
    
    selectInput("group", label = "Group in column",
                choices = unique(endPointInfo[,input$groupCol])[!is.na(unique(endPointInfo[,input$groupCol]))],
                selected = unique(endPointInfo[,10])[5],
                multiple = FALSE)
  })
  
  output$table <- DT::renderDataTable({
    statsOfSum()
  }, options = list(pageLength = 10))
  
  output$graph <- renderPlot({
    
    endpointSummary <- endpointSummary()
    
    if(nrow(endpointSummary) == 0){
      endpointSummary$site <- stationINFO$fullSiteID
    }
    
    siteLimits <- stationINFO %>%
      mutate(lakeCat = factor(Lake, 
                              levels=c("Lake Superior","Lake Michigan",
                                       "Lake Huron", "St. Lawrence River",
                                       "Detroit River and Lake St. Clair","Lake Erie","Lake Ontario"))) %>%
      arrange(lakeCat, dec.long.va) %>%
      mutate(lakeColor = c("red","black","green","brown","brown","brown","blue")[as.numeric(lakeCat)] )%>%
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
      geom_text(data=data.frame(), aes(x=c(5, 18,31,45,56),
                                       y=-.5, label=c("Superior","Michigan","Huron","Erie","Ontario")),
                colour=factor(c("red","black","green","brown","blue"),
                              levels=c("red","black","green","brown","blue")), size=3)
    print(sToxWS)
  })
  
  output$mymap <- leaflet::renderLeaflet({
    
    mapData <- left_join(stationINFO, summary, by=c("shortName"="site"))
    mapData <- mapData[!is.na(mapData$dec.lat.va),]
    mapData <- mapData[!is.na(mapData$nChem),]
    
    
    col_types <- c("darkblue","dodgerblue","green","yellow","orange","red")
    leg_vals <- c(0,1,5,10,100,1000,2500)#seq(0,2500, 500) # This is maxEAR 
    
    cols <- colorNumeric(col_types, domain = leg_vals)
    rad <- seq(5000,15000, 500)
    mapData$sizes <- rad[as.numeric(cut(mapData$nChem, c(-1:19)))]
    
    pal = colorBin(col_types, mapData$maxEAR, bins = leg_vals)
    
    leaflet(mapData) %>%
      addProviderTiles("Esri.WorldPhysical") %>%
      setView(lng = -83.5, lat = 44.5, zoom=6)%>%
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