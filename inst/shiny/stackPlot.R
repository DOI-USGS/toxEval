stackBarGroup_create <- reactive({
  catType = as.numeric(input$radioMaxGroup)
  
  chemicalSummary <- chemicalSummary()
  mean_logic <- as.logical(input$meanEAR)
  
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
                               include_legend = include_legend,
                               font_size = ifelse(catType == 2, 14, 17),
                               title = stackTitle())
  
  updateAceEditor(session, editorId = "barCode_out", value = barCode() )
  
  return(upperPlot)
})

stackTitle <- reactive({

  catType = as.numeric(input$radioMaxGroup)
  mean_logic <- as.logical(input$meanEAR)
  
  category <- c("Biological","Chemical","Chemical Class")[catType]
  
  pretty_cat <- tolower(category)
  if(pretty_cat == "biological"){
    pretty_cat <- "biological activity grouping"
  }
  title <- paste(ifelse(mean_logic,"Mean","Maximum"),"EAR",
                 "grouped by", pretty_cat)
  site <- input$sites
  browser()
  siteTable <- rawData()[["chem_site"]]
  if(site != "All"){
    title <- paste("Individual samples grouped by",pretty_cat,"
                   ",siteTable[["Fullname"]][which(siteTable$`Short Name` == site)])
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
  
  stackPlotCode <- paste0(rCodeSetup(),"
# To re-order the x-axis, 
# Convert tox_list$chem_site$`Short Name` to a factor,
# and re-order the 'levels' of that factor
stack_plot <- plot_tox_stacks(chemicalSummary, 
                  chem_site = tox_list$chem_site,
                  category = '",category,"',
                  mean_logic = ",as.logical(input$meanEAR),",
                  title = '",stackTitle(),"',
                  include_legend = ",include_legend,")
gb <- ggplot2::ggplot_build(stack_plot)
gt <- ggplot2::ggplot_gtable(gb)

gt$layout$clip[gt$layout$name=='panel-1-1'] = 'off'

grid::grid.draw(gt)")
  
  return(stackPlotCode)
  
})