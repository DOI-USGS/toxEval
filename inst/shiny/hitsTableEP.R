output$hitsTableEPs <- DT::renderDataTable({
  
  validate(
    need(!is.null(input$data), "Please select a data set")
  )
  
  chemicalSummary <- chemicalSummary()
  meanEARlogic <- as.logical(input$meanEAR)
  catType <- as.numeric(input$radioMaxGroup)
  mean_logic <- as.logical(input$meanEAR)
  hitThres <- hitThresValue()
  
  tableEPs <- table_endpoint_hits(chemicalSummary, 
                                  category = c("Biological","Chemical","Chemical Class")[catType],
                                  mean_logic = mean_logic,
                                  hit_threshold = hitThres)
  tableEPs
})

output$hitsTableEPCode <- renderPrint({
  
  catType = as.numeric(input$radioMaxGroup)
  category <- c("Biological","Chemical","Chemical Class")[catType]
  hitThres <- hitThresValue()
  
  hitsTableEPCode <- paste0(rCodeSetup(),"
table_endpoint_hits(chemicalSummary, 
                    category = '",category,"',
                    mean_logic = ",as.logical(input$meanEAR),",
                    hit_threshold = ",hitThres,")")
  
  HTML(hitsTableEPCode)
  
})