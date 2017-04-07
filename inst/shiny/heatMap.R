output$graphHeat <- renderPlot({
  
  validate(
    need(!is.null(input$data), "Please select a data set")
  )
  
  catType = as.numeric(input$radioMaxGroup)
  
  chemicalSummary <- chemicalSummary()
  rawData <- rawData()
  chem_site <- rawData$chem_site
  
  if(all(unique(chem_site$site_grouping) %in% great_lakes)){
    chem_site$site_grouping <- factor(chem_site$site_grouping,
                                      levels=great_lakes)      
  }
  
  if(all(unique(chem_site$`Short Name`) %in% sitesOrdered)){
    chem_site$`Short Name` <- factor(chem_site$`Short Name`,
                                     levels=sitesOrdered[sitesOrdered %in% unique(chem_site$`Short Name`)])
  }
  
  heatMap <- plot_tox_heatmap(chemicalSummary,
                              chem_site,
                              category = c("Biological","Chemical","Chemical Class")[catType])
  
  ggsave("heatPlot.png",heatMap,bg = "transparent")
  
  heatMap
  
})

output$graphHeat.ui <- renderUI({
  heightOfGraph <- 500
  if(as.numeric(input$radioMaxGroup) == 2){
    heightOfGraph <- 800
  }
  plotOutput("graphHeat", height = heightOfGraph)
})

output$downloadHeatPlot <- downloadHandler(
  
  filename = function() {
    "heatPlot.png"
  },
  content = function(file) {
    file.copy("heatPlot.png", file)
  }
)