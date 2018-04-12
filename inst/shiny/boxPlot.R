boxPlots_create <- reactive({
  
  catType = as.numeric(input$radioMaxGroup)
  
  plot_ND = input$plot_ND
  chemicalSummary <- chemicalSummary()
  category <- c("Biological","Chemical","Chemical Class")[catType]
  mean_logic <- as.logical(input$meanEAR)
  site <- input$sites
  
  pretty_cat <- tolower(category)
  if(pretty_cat == "biological"){
    pretty_cat <- "biological activity grouping"
  }
  
  title <- paste(ifelse(mean_logic,"Mean","Maximum"),"EAR",
                 "per site, grouped by", pretty_cat)
  if(site != "All"){
    title <- paste(title,"
                   ",rawData()[["chem_site"]][["Fullname"]])
  }

  
  bioPlot <- plot_tox_boxplots(chemicalSummary, 
                               category = category,
                               mean_logic = mean_logic,
                               plot_ND = plot_ND,
                               font_size = 18,
                               title = title)
  return(bioPlot)
  
})

boxPlot_prints <- reactive({
  
  catType = as.numeric(input$radioMaxGroup)
  
  plot_ND = input$plot_ND
  chemicalSummary <- chemicalSummary()
  category <- c("Biological","Chemical","Chemical Class")[catType]
  height <- PlotHeight()
  if(height > 750){
    text_size <- 10
  } else {
    text_size <- NA
  }
  
  bioPlot <- plot_tox_boxplots(chemicalSummary, 
                               category = category,
                               mean_logic = as.logical(input$meanEAR),
                               plot_ND = plot_ND,
                               font_size = text_size)
  return(bioPlot)
})

output$graphGroup <- renderPlot({ 

  validate(
    need(!is.null(input$data), "Please select a data set")
  )
  updateAceEditor(session, editorId = "boxCode_out", value = boxCode() )
  
  gb <- ggplot2::ggplot_build(boxPlots_create())
  gt <- ggplot2::ggplot_gtable(gb)

  gt$layout$clip[gt$layout$name=="panel"] <- "off"

  grid::grid.draw(gt)
  
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
    ggsave(file, plot = boxPlot_prints(), device = device)
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
  
  plot_ND = input$plot_ND
  
  category <- c("Biological","Chemical","Chemical Class")[catType]
  
  bioPlotCode <- paste0(rCodeSetup(),"
bio_plot <- plot_tox_boxplots(chemicalSummary, 
                  category = '",category,"',
                  mean_logic = ",as.logical(input$meanEAR),",
                  plot_ND = ",plot_ND,")
gb <- ggplot2::ggplot_build(bio_plot)
gt <- ggplot2::ggplot_gtable(gb)
gt$layout$clip[gt$layout$name=='panel'] <- 'off'
grid::grid.draw(gt)")
  
  return(bioPlotCode)
  
})

