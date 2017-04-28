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

output$endpointGraph.ui <- renderUI({
  
  height <- PlotHeight_ep()
  
  plotOutput("endpointGraph", height = height, width = 1000)
})

PlotHeight_ep = reactive({
  
  filterBy <- input$epGroup
  catType = as.numeric(input$radioMaxGroup)
  cat_col <- c("Bio_category","chnm","Class")[catType]
  
  chemicalSummary <- chemicalSummary()
  
  if(filterBy != "All"){
    chemicalSummary <- chemicalSummary %>%
      filter_(paste0(cat_col," == '", filterBy,"'"))
  }
  
  n <- 35*length(unique(chemicalSummary$endPoint))
  
  if(n < 500){
    return(500)
  } else {
    return(n)
  }
  
})

output$downloadEndpoint <- downloadHandler(
  
  filename = function() {
    "endPoint.png"
  },
  content = function(file) {
    file.copy("endPoint.png", file)
  }
)