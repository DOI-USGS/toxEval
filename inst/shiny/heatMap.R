heatMap_create <- reactive({
  plot_ND = input$plot_ND_heat
  catType = as.numeric(input$radioMaxGroup)
  
  chemicalSummary <- chemicalSummary()
  rawData <- rawData()
  chem_site <- rawData$chem_site

  if("site_grouping" %in% names(chem_site) && all(unique(chem_site$site_grouping) %in% great_lakes)){
    chem_site$site_grouping <- factor(chem_site$site_grouping,
                                      levels=great_lakes)      
  }
  
  if(all(unique(chem_site$`Short Name`) %in% sitesOrdered)){
    chem_site$`Short Name` <- factor(chem_site$`Short Name`,
                                     levels=sitesOrdered[sitesOrdered %in% unique(chem_site$`Short Name`)])
  }
  
  heatMap <- plot_tox_heatmap(chemicalSummary,
                              chem_site,
                              category = c("Biological","Chemical","Chemical Class")[catType],
                              plot_ND = plot_ND)
  return(heatMap)
})

output$graphHeat <- renderPlot({
  
  validate(
    need(!is.null(input$data), "Please select a data set")
  )
  
  print(heatMap_create())
  
})

output$graphHeat.ui <- renderUI({
  height <- PlotHeight()
  
  plotOutput("graphHeat", height = height, width="100%")
})

output$downloadHeatPlot <- downloadHandler(
  
  filename = "heatPlot.png",

  content = function(file) {
    device <- function(..., width, height) {
      grDevices::png(..., width = width, height = height,
                     res = 300, units = "in")
    }
    ggsave(file, plot = heatMap_create(), device = device)
  }
)

output$downloadHeatPlot_csv <- downloadHandler(
  
  filename = "heatPlot.csv",
  
  content = function(file) {
    write.csv(heatMap_create()[['data']], file = file, row.names = FALSE)
  }
)