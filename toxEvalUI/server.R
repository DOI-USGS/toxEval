library(DT)
library(dplyr)
library(toxEval)
library(magrittr)
library(ggplot2)

shinyServer(function(input, output) {
  
  endpointSummary <- reactive({
    endPointInfoSub <- unique(filter_(endPointInfo, 
                                      paste0(input$groupCol," == '", 
                                             ifelse(is.null(input$group),endPointInfo[1,11],input$group),
                                             "'")))
    
    endpointSummary <- chemicalSummary %>%
      filter(endPoint %in% endPointInfoSub$assay_component_endpoint_name ) %>%
      left_join(endPointInfoSub, by=c("endPoint"="assay_component_endpoint_name")) %>%
      select_("hits","EAR","endPoint","site","date",input$groupCol) %>%
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
    siteLimits <- stationINFO %>%
      mutate(lakeCat = factor(Lake, 
                              levels=c("Lake Superior","Lake Michigan",
                                       "Lake Huron", "St. Lawrence River",
                                       "Detroit River and Lake St. Clair","Lake Erie","Lake Ontario"))) %>%
      arrange(lakeCat) %>%
      mutate(lakeColor = c("red","black","green","brown","brown","brown","blue")[as.numeric(lakeCat)] )%>%
      filter(fullSiteID %in% unique(endpointSummary()$site))
    
    endPointSummBP <- endpointSummary() %>%
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

})