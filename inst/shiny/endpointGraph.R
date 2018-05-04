endpointGraph_create <- reactive({

  filterBy <- epDF[['epGroup']]
  hitThres <- hitThresValue()
  chemicalSummary <- chemicalSummary()
  catType <- as.numeric(input$radioMaxGroup)
  mean_logic <- as.logical(input$meanEAR)
  sum_logic <- as.logical(input$sumEAR)
  category <- c("Biological","Chemical","Chemical Class")[catType]

  endpointGraph <- plot_tox_endpoints(chemicalSummary, 
                                      category = category,
                                      mean_logic = mean_logic,
                                      sum_logic = sum_logic,
                                      hit_threshold = hitThresValue(),
                                      filterBy = filterBy,
                                      title = epTitle(),
                                      font_size = 18) 
  
  updateAceEditor(session, editorId = "epGraph_out", value = epGraphCode() )
  return(endpointGraph)
})

output$endpointGraph <- renderPlot({ 
  
  validate(
    need(!is.null(input$data), "Please select a data set")
  )
  
  gb <- ggplot2::ggplot_build(endpointGraph_create())
  gt <- ggplot2::ggplot_gtable(gb)
  gt$layout$clip[gt$layout$name=="panel"] <- "off"
  grid::grid.draw(gt)

})

output$endpointGraph.ui <- renderUI({
  
  height <- PlotHeight_ep()
  
  plotOutput("endpointGraph", height = height, width = "100%")
})

PlotHeight_ep = reactive({
  
  filterBy <- epDF[['epGroup']]
  catType = as.numeric(input$radioMaxGroup)
  cat_col <- c("Bio_category","chnm","Class")[catType]
  
  chemicalSummary <- chemicalSummary()
  
  if(filterBy != "All"){
    chemicalSummary <- chemicalSummary %>%
      filter_(paste0(cat_col," == '", filterBy,"'"))
  }
  
  n <- 35*length(unique(chemicalSummary$endPoint))
  
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

epTitle <- reactive({
  catType <- as.numeric(input$radioMaxGroup)
  mean_logic <- as.logical(input$meanEAR)
  sum_logic <- as.logical(input$sumEAR)
  site <- input$sites
  siteTable <- rawData()[["chem_site"]]
  category <- c("Biological","Chemical","Chemical Class")[catType]
  filterBy <- epDF[['epGroup']]
 
  if(site == "All"){
    pretty_cat <- switch(category, 
                         "Chemical" = paste("for",filterBy),
                         "Biological" = paste("for chemicals within the",filterBy,"class"),
                         "Chemical Class" = paste("for chemicals within the",filterBy,"group")
    )
    pretty_cat <- paste(pretty_cat,"
")
    title <- paste("Summing EARs",pretty_cat,"for a given sample,")
    if(mean_logic){
      title <- paste(title, "taking the mean of each site")
    } else {
      title <- paste(title, "taking the max of each site")
    }
  } else {
    pretty_cat <- switch(category, 
                         "Chemical" = filterBy,
                         "Biological" = paste("chemicals within",filterBy),
                         "Chemical Class" = paste("chemicals within",filterBy)
    )
    title <- paste0("EAR per ",category)
      
    title <- paste(title,"
", siteTable[["Fullname"]][which(siteTable$`Short Name` == site)])
  }

  return(title)
})

epGraphCode <- reactive({
  
  catType = as.numeric(input$radioMaxGroup)
  category <- c("Biological","Chemical","Chemical Class")[catType]
  hitThres <- hitThresValue()
  filterBy <- epDF[['epGroup']]
  mean_logic <- as.logical(input$meanEAR)
  sum_logic <- as.logical(input$sumEAR)
  if(sum_logic){
    epGraphCode <- paste0(rCodeSetup(),"
ep_plot <- plot_tox_endpoints(chemicalSummary, 
                    category = '",category,"',
                    mean_logic = ",mean_logic,",
                    hit_threshold = ",hitThres,",
                    title = '",epTitle(),"'
                    filterBy = '",filterBy,"')
ep_plot")
  } else {
    epGraphCode <- paste0(rCodeSetup(),"
ep_plot <- plot_tox_endpoints(chemicalSummary, 
                        category = '",category,"',
                        mean_logic = ",mean_logic,",
                        sum_logic = FALSE,
                        hit_threshold = ",hitThres,",
                        title = '",epTitle(),"'
                        filterBy = '",filterBy,"')
ep_plot")    
  }
  
  return(epGraphCode)
  
})