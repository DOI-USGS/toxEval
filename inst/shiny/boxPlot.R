output$graphGroup <- renderPlot({ 
  
  catType = as.numeric(input$radioMaxGroup)

  plot_ND = input$plot_ND
  
  chemicalSummary <- chemicalSummary()
  
  validate(
    need(!is.null(input$data), "Please select a data set")
  )

  category <- c("Biological","Chemical","Chemical Class")[catType]

  bioPlot <- plot_tox_boxplots(chemicalSummary, 
                               category = category,
                               mean_logic = as.logical(input$meanEAR),
                               plot_ND = plot_ND)
  
  if(catType == 2){
    ggsave("boxPlot.png",
           bioPlot,
           bg = "transparent")
  } else {
    ggsave("boxPlot.png",
           bioPlot,
           bg = "transparent")
  }
  
  bioPlot
  
})

PlotHeight = reactive({
  
  catType = as.numeric(input$radioMaxGroup)
  cat_col <- c("Bio_category","chnm","Class")[catType]
  chemicalSummary <- chemicalSummary()
  
  n <- 35*length(unique(chemicalSummary[[cat_col]]))

  if(n < 500){
    return(500)
  } else {
    return(n)
  }

})

output$graphGroup.ui <- renderUI({
  
  height <- PlotHeight()

  plotOutput("graphGroup", height = height, width = 1000)
})

output$downloadBoxPlot <- downloadHandler(
  
  filename = function() {
    "boxPlot.png"
  },
  content = function(file) {
    file.copy("boxPlot.png", file)
  }
)