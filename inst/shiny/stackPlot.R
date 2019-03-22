stackBarGroup_create <- reactive({
  catType = as.numeric(input$radioMaxGroup)

  chemical_summary <- chemical_summary()
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

  top_num <- ifelse(catType == 2, 5, NA)

  category <- c("Biological","Chemical","Chemical Class")[catType]

  upperPlot <- plot_tox_stacks(chemical_summary, 
                               chem_site, 
                               category = category,
                               mean_logic = mean_logic,
                               sum_logic = sum_logic,
                               include_legend = TRUE,
                               font_size = ifelse(catType == 2, 14, 17),
                               title = genericTitle(),
                               top_num = top_num)
  
  shinyAce::updateAceEditor(session, editorId = "barCode_out", value = barCode() )
  
  return(upperPlot)
})

output$stackBarGroup <- renderPlot({
  
  validate(
  need(!is.null(rawData_data$data), "Please select a data set")
  )
  
  stackBarGroup_create()

})

output$downloadStackPlot <- downloadHandler(
  
  filename = "stackPlot.png",
  
  content = function(file) {
    ggplot2::ggsave(file, plot = stackBarGroup_create(),  
                    device = "png", width = 11,
                    height = 9)
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
  top_num <- ifelse(catType == 2, 5, NA)
  mean_logic <- as.logical(input$meanEAR)
  sum_logic <- as.logical(input$sumEAR)

  stackPlotCode <- paste0(rCodeSetup(),"
# To re-order the x-axis, 
# Convert tox_list$chem_site$`Short Name` to a factor,
# and re-order the 'levels' of that factor
stack_plot <- plot_tox_stacks(chemical_summary, 
                  chem_site = tox_list$chem_site,
                  category = '",category,"',
                  mean_logic = ",mean_logic,",
                  title = '",genericTitle(),"',
                  include_legend = TRUE,
                  top_num = ",top_num)

  if(!sum_logic){
    stackPlotCode <- paste0(stackPlotCode,"
                            sum_logic = FALSE)")
  } else {
    stackPlotCode <- paste0(stackPlotCode,")")
  }
    
  stackPlotCode <- paste0(stackPlotCode,"
stack_plot
# To save:
# Fiddle with height and width (in inches) for best results:
# Change file name extension to save as png.
# ggplot2::ggsave(stack_plot, file='stack.pdf',
#                        height = 5,
#                        width = 7)")
  
  return(stackPlotCode)
  
})