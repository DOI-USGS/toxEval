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

output$tableSummCode <- renderPrint({
  
  catType = as.numeric(input$radioMaxGroup)
  category <- c("Biological","Chemical","Chemical Class")[catType]
  hitThres <- hitThresValue()
  
  tableSummCode <- paste0(rCodeSetup(),"
table_tox_rank(chemicalSummary, 
                  category = '",category,"',
                  mean_logic = ",as.logical(input$meanEAR),",
                  hit_threshold = ",hitThres,")")
  
  HTML(tableSummCode)
  
})