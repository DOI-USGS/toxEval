output$endpointGraph <- renderPlot({ 
  
  filterBy <- input$epGroup
  meanEARlogic <- as.logical(input$meanEAR)
  hitThres <- hitThresValue()
  chemicalSummary <- chemicalSummary()
  catType = as.numeric(input$radioMaxGroup)
  
  endpointGraph <- plot_tox_endpoints(chemicalSummary, 
                                      category = c("Biological","Chemical","Chemical Class")[catType],
                                      mean_logic = as.logical(input$meanEAR),
                                      hit_threshold = hitThresValue(),
                                      filterBy = filterBy) 
  
  ggsave("endPoint.png",
         endpointGraph,
         bg = "transparent", 
         height = 12, width = 9)
  
  endpointGraph
})

output$downloadEndpoint <- downloadHandler(
  
  filename = function() {
    "endPoint.png"
  },
  content = function(file) {
    file.copy("endPoint.png", file)
  }
)