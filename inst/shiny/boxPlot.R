boxPlots_create <- reactive({
  
  catType = as.numeric(input$radioMaxGroup)
  
  plot_ND = input$plot_ND
  
  chemicalSummary <- chemicalSummary()
  category <- c("Biological","Chemical","Chemical Class")[catType]
  
  bioPlot <- plot_tox_boxplots(chemicalSummary, 
                               category = category,
                               mean_logic = as.logical(input$meanEAR),
                               plot_ND = plot_ND)
  return(bioPlot)
  
})

output$graphGroup <- renderPlot({ 

  validate(
    need(!is.null(input$data), "Please select a data set")
  )

  print(boxPlots_create())
  
})

PlotHeight = reactive({
  
  catType = as.numeric(input$radioMaxGroup)
  cat_col <- c("Bio_category","chnm","Class")[catType]
  chemicalSummary <- chemicalSummary()
  
  n <- 35*length(unique(chemicalSummary[[cat_col]]))

  if(n < 750){
    return(750)
  } else {
    return(n)
  }

})

output$graphGroup.ui <- renderUI({
  
  height <- PlotHeight()

  plotOutput("graphGroup", height = height, width="100%")
})

output$downloadBoxPlot <- downloadHandler(
  
  filename = "boxPlot.png",
  
  content = function(file) {
    device <- function(..., width, height) {
      grDevices::png(..., width = width, height = height,
                     res = 300, units = "in")
    }
    ggsave(file, plot = boxPlots_create(), device = device)
  }
    
)