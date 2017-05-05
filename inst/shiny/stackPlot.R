stackBarGroup_create <- reactive({
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
  return(upperPlot)
})

output$stackBarGroup <- renderPlot({
  
 validate(
  need(!is.null(input$data), "Please select a data set")
 )
  
 print(stackBarGroup_create())
 
})

output$downloadStackPlot <- downloadHandler(
  
  filename = "stackPlot.png",
  
  content = function(file) {
    device <- function(..., width, height) {
      grDevices::png(..., width = width, height = height,
                     res = 300, units = "in")
    }
    ggsave(file, plot = stackBarGroup_create(), device = device)
  }
)