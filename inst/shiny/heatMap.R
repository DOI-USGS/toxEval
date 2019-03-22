heatMap_create <- reactive({
  plot_ND = input$plot_ND_heat
  catType = as.numeric(input$radioMaxGroup)
  
  chemical_summary <- chemical_summary()
  rawData <- rawData()
  chem_site <- rawData$chem_site
  mean_logic <- as.logical(input$meanEAR)
  sum_logic <- as.logical(input$sumEAR)

  if("site_grouping" %in% names(chem_site) && all(unique(chem_site$site_grouping) %in% great_lakes)){
    chem_site$site_grouping <- factor(chem_site$site_grouping,
                                      levels=great_lakes)      
  }
  
  if(all(unique(chem_site$`Short Name`) %in% sitesOrdered)){
    chem_site$`Short Name` <- factor(chem_site$`Short Name`,
                                     levels=sitesOrdered[sitesOrdered %in% unique(chem_site$`Short Name`)])
  }
  category <-  c("Biological","Chemical","Chemical Class")[catType]
  
  heatMap <- plot_tox_heatmap(chemical_summary,
                              chem_site,
                              category = category,
                              plot_ND = plot_ND,
                              mean_logic = mean_logic,
                              sum_logic = sum_logic,
                              font_size = ifelse(catType == 2, 14, 17),
                              title = genericTitle())
  
  shinyAce::updateAceEditor(session, editorId = "heat_out", value = heatCode() )
  
  return(heatMap)
})

output$graphHeat <- renderPlot({
  
  validate(
    need(!is.null(rawData_data$data), "Please select a data set")
  )
  
  print(heatMap_create())
  
})

output$graphHeat.ui <- renderUI({
  height <- PlotHeight()
  
  shinycssloaders::withSpinner(plotOutput("graphHeat", height = height, width="100%"))
})

output$downloadHeatPlot <- downloadHandler(
  
  filename = "heatPlot.png",

  content = function(file) {
    ggplot2::ggsave(file, plot = heatMap_create(), 
                    device = "png", width = 11,
                    height = PlotHeight()/200)
  }
)

output$downloadHeatPlot_csv <- downloadHandler(
  
  filename = "heatPlot.csv",
  
  content = function(file) {
    write.csv(heatMap_create()[['data']], file = file, row.names = FALSE)
  }
)

heatCode <- reactive({
  
  catType = as.numeric(input$radioMaxGroup)
  category <- c("Biological","Chemical","Chemical Class")[catType]
  plot_ND = input$plot_ND_heat
  mean_logic <- as.logical(input$meanEAR)
  sum_logic <- as.logical(input$sumEAR)
  
  heatCode <- paste0(rCodeSetup(),"
# To re-order the x-axis, 
# Convert tox_list$chem_site$`Short Name` to a factor,
# and re-order the 'levels' of that factor
heat_map <- plot_tox_heatmap(chemical_summary,
            chem_site = tox_list$chem_site,
            category = '",category,"',
            mean_logic = ",mean_logic,",
            title = '",genericTitle(),"',
            plot_ND = ",plot_ND)
  if(!sum_logic){
    heatCode <- paste0(heatCode,",
            sum_logic = FALSE)")
  } else {
    heatCode <- paste0(heatCode,")")
  }
  
  heatCode <- paste0(heatCode,"
heat_map
# To save:
# Fiddle with height and width (in inches) for best results:
                     # Change file name extension to save as png.
                     # ggplot2::ggsave(heat_map, file='heat.pdf',
                     #                        height = 11,
                     #                        width = 7)")
  
  return(heatCode)
  
})