library(DT)
library(dplyr)
library(toxEval)

shinyServer(function(input, output) {
  
  output$groupControl <- renderUI({
    
    selectInput("group", label = "Group in column",
                choices = unique(endPointInfo[,input$groupCol])[!is.na(unique(endPointInfo[,input$groupCol]))],
                selected = unique(endPointInfo[,10])[5],
                multiple = FALSE)
  })
  
  output$table <- renderDataTable({

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
    
    statsOfSum <- endpointSummary%>%
      group_by(site) %>%
      summarise(meanEAR = mean(sumEAR),
                 maxEAR = max(sumEAR),
                 sumHits = sum(nHits),
                 nSamples = n()) %>%
      data.frame()%>%
      mutate(site = siteKey[site])
    
      statsOfSum
  }, options = list(pageLength = 10)) 
})