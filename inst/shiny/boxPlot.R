output$graphGroup <- renderPlot({ 
  
  catType = as.numeric(input$radioMaxGroup)
  
  chemicalSummary <- chemicalSummary()
  
  validate(
    need(!is.null(input$data), "Please select a data set")
  )
  
  bioPlot <- plot_tox_boxplots(chemicalSummary, 
                               category = c("Biological","Chemical","Chemical Class")[catType],
                               mean_logic = as.logical(input$meanEAR))
  
  if(catType == 2){
    ggsave("boxPlot.png",
           bioPlot,
           bg = "transparent", 
           height = 12, width = 9)
  } else {
    ggsave("boxPlot.png",
           bioPlot,
           bg = "transparent", 
           height = 6, width = 8)
  }
  
  bioPlot
  
})

output$graphGroup.ui <- renderUI({
  heightOfGraph <- 500
  if(as.numeric(input$radioMaxGroup) == 2){
    heightOfGraph <- 800
  }
  plotOutput("graphGroup", height = heightOfGraph)
})

output$downloadBoxPlot <- downloadHandler(
  
  filename = function() {
    "boxPlot.png"
  },
  content = function(file) {
    file.copy("boxPlot.png", file)
  }
)