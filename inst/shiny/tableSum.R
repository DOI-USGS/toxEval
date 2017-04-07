output$tableSumm <- DT::renderDataTable({
  
  validate(
    need(!is.null(input$data), "Please select a data set")
  )
  
  catType = as.numeric(input$radioMaxGroup)
  
  chemicalSummary <- chemicalSummary()
  hitThres <- hitThresValue()
  mean_logic <- as.logical(input$meanEAR)
  
  tableGroup <- table_tox_rank(chemicalSummary, 
                               category = c("Biological","Chemical","Chemical Class")[catType],
                               mean_logic = mean_logic,
                               hit_threshold = hitThres)
  tableGroup
  
})