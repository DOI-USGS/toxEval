library(DT)
library(dplyr)
library(toxEval)

shinyServer(function(input, output) {
  
  output$groupControl <- renderUI({
    
    selectInput("group", label = "Group in column",
                choices = unique(endPointInfo[,input$groupCol])[!is.na(unique(endPointInfo[,input$groupCol]))],
                multiple = FALSE)
  })
  
  output$table <- renderDataTable({
    
    endPointInfoSub <- unique(filter_(endPointInfo, paste0(input$groupCol," == '", input$group, "'")))
    
    endpointSummary <- chemicalSummary %>%
      filter(endPoint %in% endPointInfoSub$assay_component_endpoint_name ) %>%
      left_join(endPointInfoSub, by=c("endPoint"="assay_component_endpoint_name")) %>%
      select_("hits","EAR","endPoint","site","date",input$groupCol) %>%
      group_by(site,date) %>%
      summarise(sumEAR = sum(EAR),
                nHits = sum(hits)) 
    
    statsOfSum <- endpointSummary%>%
      group_by_(chemicalSummary.site) %>%
      summarise_(meanEAR = paste0("mean(",funcToSummerise,")"),
                 # medianEAR = paste0("median(",funcToSummerise,")"),
                 maxEAR = paste0("max(",funcToSummerise,")"),
                 sumHits = "sum(nHits)",
                 nSamples = "n()") %>%
      data.frame()%>%
      mutate(site = siteKey[site])
    
      statsOfSum
  })
})