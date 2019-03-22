boxPlots_create <- reactive({
  
  catType = as.numeric(input$radioMaxGroup)
  include_thresh <- as.logical(input$plot_thres_box)
  hitThres <- ifelse(include_thresh, hitThresValue(),NA)
  plot_ND = input$plot_ND
  chemical_summary <- chemical_summary()
  category <- c("Biological","Chemical","Chemical Class")[catType]

  mean_logic <- as.logical(input$meanEAR)
  sum_logic <- as.logical(input$sumEAR)

  bioPlot <- plot_tox_boxplots(chemical_summary, 
                               category = category,
                               mean_logic = mean_logic,
                               sum_logic = sum_logic,
                               plot_ND = plot_ND,
                               hit_threshold = hitThres,
                               font_size = 18,
                               title = genericTitle())
  return(bioPlot)
  
})

boxPlot_prints <- reactive({
  
  catType = as.numeric(input$radioMaxGroup)
  
  plot_ND = input$plot_ND
  chemical_summary <- chemical_summary()
  category <- c("Biological","Chemical","Chemical Class")[catType]
  height <- PlotHeight()
  if(height > 750){
    text_size <- 10
  } else {
    text_size <- NA
  }
  
  bioPlot <- plot_tox_boxplots(chemical_summary, 
                               category = category,
                               mean_logic = as.logical(input$meanEAR),
                               sum_logic = as.logical(input$sumEAR),
                               plot_ND = plot_ND,
                               font_size = text_size)
  return(bioPlot)
})

output$graphGroup <- renderPlot({ 
  
  validate(
    need(!is.null(rawData_data$data), "Please select a data set")
  )
  shinyAce::updateAceEditor(session, editorId = "boxCode_out", value = boxCode() )
  
  boxPlots_create()
  
}) 

PlotHeight = reactive({
  
  catType = as.numeric(input$radioMaxGroup)
  cat_col <- c("Bio_category","chnm","Class")[catType]
  chemical_summary <- chemical_summary()
  
  n <- 35*length(unique(chemical_summary[[cat_col]]))

  if(n < 750){
    return(750)
  } else {
    return(n)
  }

})

output$graphGroup.ui <- renderUI({
  
  height <- PlotHeight()

  shinycssloaders::withSpinner(plotOutput("graphGroup", height = height, width="100%"))
})

output$downloadBoxPlot <- downloadHandler(

  filename = "boxPlot.png",
  content = function(file) {
    ggplot2::ggsave(file, plot = boxPlot_prints(), 
                    device = "png", width = 11,
                    height = PlotHeight()/200)
  }
)

output$downloadBoxPlot_csv <- downloadHandler(
  
  filename = "boxPlot.csv",
  
  content = function(file) {
    write.csv(boxPlots_create()[['data']], file, row.names = FALSE)
  }
)

boxCode <- reactive({
  
  catType = as.numeric(input$radioMaxGroup)
  mean_logic <- as.logical(input$meanEAR)
  sum_logic <- as.logical(input$sumEAR)
  plot_ND = input$plot_ND
  
  include_thresh <- as.logical(input$plot_thres_box)
  hitThres <- ifelse(include_thresh, hitThresValue(),NA)
  category <- c("Biological","Chemical","Chemical Class")[catType]
  
bioPlotCode <- paste0(rCodeSetup(),"
bio_plot <- plot_tox_boxplots(chemical_summary, 
                        category = '",category,"',
                        mean_logic = ",mean_logic,",
                        hit_threshold = ",hitThres,",
                        title = '",genericTitle(),"',
                        plot_ND = ",plot_ND)
  if(sum_logic){
    bioPlotCode <- paste0(bioPlotCode,")")
  } else {
    bioPlotCode <- paste0(bioPlotCode,"
                          sum_logic = FALSE)")

  }

  bioPlotCode <- paste0(bioPlotCode,"
bio_plot
# To save:
# Fiddle with height and width (in inches) for best results:
# Change file name extension to save as png.
# ggplot2::ggsave(bio_plot, file='boxplot.pdf',
#                        height = 9,
#                        width = 11)")
  
  return(bioPlotCode)
  
})

