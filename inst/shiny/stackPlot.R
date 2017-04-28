output$downloadStackPlot <- downloadHandler(
  
  filename = function() {
    "stackPlot.png"
  },
  content = function(file) {
    file.copy("stackPlot.png", file)
  }
)

output$stackBarGroup <- renderPlot({
  
  validate(
    need(!is.null(input$data), "Please select a data set")
  )
  
  catType = as.numeric(input$radioMaxGroup)
  
  chemicalSummary <- chemicalSummary()
  mean_logic <- as.logical(input$meanEAR)
  
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
  
  upperPlot <- plot_tox_stacks(chemicalSummary, 
                               chem_site, 
                               category = c("Biological","Chemical","Chemical Class")[catType],
                               mean_logic = mean_logic)
  
  ggsave("stackPlot.png",upperPlot,bg = "transparent")
  
  upperPlot
}, height = 600, width = 1000)