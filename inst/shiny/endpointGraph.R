endpointGraph_create <- reactive({

  filterBy <- epDF[['epGroup']]
  include_thresh <- as.logical(input$plot_thres_ep)
  hitThres <- ifelse(include_thresh, hitThresValue(),NA)
  chemical_summary <- chemical_summary()
  catType <- as.numeric(input$radioMaxGroup)
  top_num <- as.numeric(input$topNum)
  mean_logic <- as.logical(input$meanEAR)
  sum_logic <- as.logical(input$sumEAR)
  category <- c("Biological","Chemical","Chemical Class")[catType]

  validate(
    need(top_num <= 50 , "Shiny app cannot display more than 50 endpoints")
  )
  
  endpointGraph <- plot_tox_endpoints(chemical_summary, 
                                      category = category,
                                      mean_logic = mean_logic,
                                      sum_logic = sum_logic,
                                      hit_threshold = hitThres,
                                      filterBy = filterBy,
                                      title = genericTitle(),
                                      font_size = 18,
                                      top_num = top_num) 
  
  shinyAce::updateAceEditor(session, editorId = "epGraph_out", value = epGraphCode() )
  return(endpointGraph)
})

output$endpointGraph <- renderPlot({ 
  
  validate(
    need(!is.null(rawData_data$data), "Please select a data set")
  )
  
  endpointGraph_create()

})

output$endpointGraph.ui <- renderUI({
  
  height <- PlotHeight_ep()
  
  plotOutput("endpointGraph", height = height, width = "100%")
})

PlotHeight_ep = reactive({
  
  filterBy <- epDF[['epGroup']]
  catType = as.numeric(input$radioMaxGroup)
  cat_col <- c("Bio_category","chnm","Class")[catType]
  
  top_num <- as.numeric(input$topNum)
  
  chemical_summary <- chemical_summary()
  
  if(filterBy != "All"){
    chemical_summary <- chemical_summary[chemical_summary[cat_col] == filterBy,]
  }
  
  if(is.na(top_num)){
    n <- 35*length(unique(chemical_summary$endPoint))
  } else {
    n <- 35*top_num
  }
  
  if(n < 500){
    return(500)
  } else {
    return(n)
  }
  
})

output$downloadEndpoint <- downloadHandler(
  
  filename = "endPoint.png",
  
  content = function(file) {
    ggplot2::ggsave(file, plot = endpointGraph_create(), 
                    device = "png", width = 11,
                    height = PlotHeight_ep()/300)
  }
)

output$downloadEndpoint_csv <- downloadHandler(
  
  filename = "endPoint.csv",
  
  content = function(file) {
    write.csv(endpointGraph_create()[['data']], file, row.names = FALSE)
  }
)

epGraphCode <- reactive({
  
  catType = as.numeric(input$radioMaxGroup)
  category <- c("Biological","Chemical","Chemical Class")[catType]
  include_thresh <- as.logical(input$plot_thres_ep)
  hitThres <- ifelse(include_thresh, hitThresValue(),NA)
  filterBy <- epDF[['epGroup']]
  mean_logic <- as.logical(input$meanEAR)
  sum_logic <- as.logical(input$sumEAR)
  top_num <- as.numeric(input$topNum)
  
  epGraphCode <- paste0(rCodeSetup(),"
ep_plot <- plot_tox_endpoints(chemical_summary, 
                        category = '",category,"',
                        mean_logic = ",mean_logic,",
                        hit_threshold = ",hitThres,",
                        title = '",genericTitle(),"',
                        top_num = ",top_num,",
                        filterBy = '",filterBy,"'")
  
  if(sum_logic){
    epGraphCode <- paste0(epGraphCode,")")

  } else {
    epGraphCode <- paste0(epGraphCode,",
                        sum_logic = FALSE)")
  }
  epGraphCode <- paste0(epGraphCode,"  
ep_plot
# To save:
# Fiddle with height and width (in inches) for best results:
# Change file name extension to save as png.
# ggplot2::ggsave(ep_plot, file='ep_box.pdf',
#                        height = 7,
#                        width = 7)")
  
  return(epGraphCode)
  
})