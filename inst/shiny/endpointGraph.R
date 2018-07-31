endpointGraph_create <- reactive({

  filterBy <- epDF[['epGroup']]
  include_thresh <- as.logical(input$plot_thres_ep)
  hitThres <- ifelse(include_thresh, hitThresValue(),NA)
  chemical_summary <- chemical_summary()
  catType <- as.numeric(input$radioMaxGroup)
  mean_logic <- as.logical(input$meanEAR)
  sum_logic <- as.logical(input$sumEAR)
  category <- c("Biological","Chemical","Chemical Class")[catType]

  endpointGraph <- plot_tox_endpoints(chemical_summary, 
                                      category = category,
                                      mean_logic = mean_logic,
                                      sum_logic = sum_logic,
                                      hit_threshold = hitThres,
                                      filterBy = filterBy,
                                      title = genericTitle(),
                                      font_size = 18) 
  
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
  
  shinycssloaders::withSpinner(plotOutput("endpointGraph", height = height, width = "100%"))
})

PlotHeight_ep = reactive({
  
  filterBy <- epDF[['epGroup']]
  catType = as.numeric(input$radioMaxGroup)
  cat_col <- c("Bio_category","chnm","Class")[catType]
  
  chemical_summary <- chemical_summary()
  
  if(filterBy != "All"){
    chemical_summary <- chemical_summary %>%
      dplyr::filter_(paste0(cat_col," == '", filterBy,"'"))
  }
  
  n <- 35*length(unique(chemical_summary$endPoint))
  
  if(n < 500){
    return(500)
  } else {
    return(n)
  }
  
})

output$downloadEndpoint <- downloadHandler(
  
  filename = "endPoint.png",
  
  content = function(file) {
    device <- function(..., width, height) {
      grDevices::png(..., width = width, height = height,
                     res = 300, units = "in")
    }
    ggsave(file, plot = endpointGraph_create(), device = device)
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
  if(sum_logic){
    epGraphCode <- paste0(rCodeSetup(),"
ep_plot <- plot_tox_endpoints(chemical_summary, 
                    category = '",category,"',
                    mean_logic = ",mean_logic,",
                    hit_threshold = ",hitThres,",
                    title = '",genericTitle(),"'
                    filterBy = '",filterBy,"')
ep_plot")
  } else {
    epGraphCode <- paste0(rCodeSetup(),"
ep_plot <- plot_tox_endpoints(chemical_summary, 
                        category = '",category,"',
                        mean_logic = ",mean_logic,",
                        sum_logic = FALSE,
                        hit_threshold = ",hitThres,",
                        title = '",genericTitle(),"'
                        filterBy = '",filterBy,"')
ep_plot")    
  }
  
  return(epGraphCode)
  
})