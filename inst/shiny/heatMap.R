heatMap_create <- reactive({
  plot_ND = input$plot_ND_heat
  catType = as.numeric(input$radioMaxGroup)
  
  chemicalSummary <- chemicalSummary()
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
  
  heatMap <- plot_tox_heatmap(chemicalSummary,
                              chem_site,
                              category = category,
                              plot_ND = plot_ND,
                              mean_logic = mean_logic,
                              sum_logic = sum_logic,
                              font_size = ifelse(catType == 2, 14, 17),
                              title = heatTitle())
  
  updateAceEditor(session, editorId = "heat_out", value = heatCode() )
  
  return(heatMap)
})

heatTitle <- reactive({
  catType = as.numeric(input$radioMaxGroup)
  category <- c("Biological","Chemical","Chemical Class")[catType]
  
  mean_logic <- as.logical(input$meanEAR)
  site <- input$sites
  siteTable <- rawData()[["chem_site"]]

  if(site == "All"){
    pretty_cat <- switch(category, 
                         "Chemical" = "for all chemicals",
                         "Biological" = "for chemicals within a grouping",
                         "Chemical Class" = "for chemicals within a class"
    )
    title <- paste("Summing EARs",pretty_cat, "for a given sample,")
    if (mean_logic){
      title <- paste(title,"taking the mean of each site")
    } else {
      title <- paste(title,"taking the max of each site")
    }
  } else {
      pretty_cat <- switch(category, 
                           "Chemical" = "chemical",
                           "Biological" = "grouping",
                           "Chemical Class" = "chemical class"
      )
      word <- ifelse(mean_logic,"Mean","Maximum")
      title <- paste(word,"EAR per",pretty_cat)
      
      title <- paste(title,"
", siteTable[["Fullname"]][which(siteTable$`Short Name` == site)])
  }
  return(title)
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

heatCode <- reactive({
  
  catType = as.numeric(input$radioMaxGroup)
  category <- c("Biological","Chemical","Chemical Class")[catType]
  plot_ND = input$plot_ND_heat
  mean_logic <- as.logical(input$meanEAR)
  sum_logic <- as.logical(input$sumEAR)
  if(sum_logic){
    heatCode <- paste0(rCodeSetup(),"
# To re-order the x-axis, 
# Convert tox_list$chem_site$`Short Name` to a factor,
# and re-order the 'levels' of that factor
plot_tox_heatmap(chemicalSummary,
                 chem_site = tox_list$chem_site,
                 category = '",category,"',
                 mean_logic = ",mean_logic,",
                 title = '",heatTitle(),"',
                 plot_ND = ",plot_ND,")")
  } else {
    heatCode <- paste0(rCodeSetup(),"
# To re-order the x-axis, 
# Convert tox_list$chem_site$`Short Name` to a factor,
# and re-order the 'levels' of that factor
plot_tox_heatmap(chemicalSummary,
                 chem_site = tox_list$chem_site,
                 category = '",category,"',
                 mean_logic = ",mean_logic,",
                 sum_logic = FALSE,
                 title = '",heatTitle(),"',
                 plot_ND = ",plot_ND,")")    
  }
  
  HTML(heatCode)
  
})