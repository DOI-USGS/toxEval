stackBarGroup_create <- reactive({
  catType = as.numeric(input$radioMaxGroup)

  chemicalSummary <- chemicalSummary()
  mean_logic <- as.logical(input$meanEAR)
  sum_logic <- as.logical(input$sumEAR)
  
  rawData <- rawData()
  chem_site <- rawData$chem_site

  text_size <- input$text_size1
  
  if("site_grouping" %in% names(chem_site) && all(unique(chem_site$site_grouping) %in% great_lakes)){
    chem_site$site_grouping <- factor(chem_site$site_grouping,
                                      levels=great_lakes)      
  }
  
  if(all(unique(chem_site$`Short Name`) %in% sitesOrdered)){
    chem_site$`Short Name` <- factor(chem_site$`Short Name`,
                                     levels=sitesOrdered[sitesOrdered %in% unique(chem_site$`Short Name`)])
  }

  include_legend <- !(catType == 2)

  category <- c("Biological","Chemical","Chemical Class")[catType]

  upperPlot <- plot_tox_stacks(chemicalSummary, 
                               chem_site, 
                               category = category,
                               mean_logic = mean_logic,
                               sum_logic = sum_logic,
                               include_legend = include_legend,
                               font_size = ifelse(catType == 2, 14, 17),
                               title = stackTitle())
  
  updateAceEditor(session, editorId = "barCode_out", value = barCode() )
  
  return(upperPlot)
})

stackTitle <- reactive({

  catType = as.numeric(input$radioMaxGroup)
  
  mean_logic <- as.logical(input$meanEAR)
  sum_logic <- as.logical(input$sumEAR)
  
  category <- c("Biological","Chemical","Chemical Class")[catType]
  
  pretty_cat <- tolower(category)
  
  if(pretty_cat == "biological"){
    pretty_cat <- "grouping"
  }
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
    } else{
      title <- paste(title,"taking the max of each site")
    }
  } else {
    pretty_cat <- switch(category, 
                         "Chemical" = "Chemical",
                         "Biological" = "Biological Activity Grouping",
                         "Chemical Class" = "Chemical Class"
    )
    word <- ifelse(mean_logic,"Mean","Maximum")
    title <- paste(word,"EAR per",pretty_cat)
    
    title <- paste(title,"
", siteTable[["Fullname"]][which(siteTable$`Short Name` == site)])
  }
  return(title)
})

output$stackBarGroup <- renderPlot({
  
  validate(
  need(!is.null(input$data), "Please select a data set")
  )
  
  gb <- ggplot2::ggplot_build(stackBarGroup_create())
  gt <- ggplot2::ggplot_gtable(gb)
  
  gt$layout$clip[gt$layout$name=="panel-1-1"] = "off"
  
  grid::grid.draw(gt)

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

output$downloadStackPlot_csv <- downloadHandler(
  
  filename = "stackPlot.csv",
  
  content = function(file) {
    write.csv(stackBarGroup_create()[['data']], file, row.names = FALSE)
  }
)

barCode <- reactive({
  
  catType = as.numeric(input$radioMaxGroup)
  category <- c("Biological","Chemical","Chemical Class")[catType]
  include_legend <- !(catType == 2)
  mean_logic <- as.logical(input$meanEAR)
  sum_logic <- as.logical(input$sumEAR)
  
  if(!sum_logic){
    stackPlotCode <- paste0(rCodeSetup(),"
# To re-order the x-axis, 
# Convert tox_list$chem_site$`Short Name` to a factor,
# and re-order the 'levels' of that factor
stack_plot <- plot_tox_stacks(chemicalSummary, 
                  chem_site = tox_list$chem_site,
                  category = '",category,"',
                  mean_logic = ",mean_logic,",
                  sum_logic = FALSE,
                  title = '",stackTitle(),"',
                  include_legend = ",include_legend,")
stack_plot")    
  } else {
  stackPlotCode <- paste0(rCodeSetup(),"
# To re-order the x-axis, 
# Convert tox_list$chem_site$`Short Name` to a factor,
# and re-order the 'levels' of that factor
stack_plot <- plot_tox_stacks(chemicalSummary, 
                  chem_site = tox_list$chem_site,
                  category = '",category,"',
                  mean_logic = ",mean_logic,",
                  title = '",stackTitle(),"',
                  include_legend = ",include_legend,")
stack_plot")    
  }

  
  return(stackPlotCode)
  
})