tableSummData <- reactive({
  catType = as.numeric(input$radioMaxGroup)
  
  chemicalSummary <- chemicalSummary()
  hitThres <- hitThresValue()
  mean_logic <- as.logical(input$meanEAR)
  
  tableGroup <- stats_of_groupings(chemicalSummary, 
                               category = c("Biological","Chemical","Chemical Class")[catType],
                               mean_logic = mean_logic,
                               hit_threshold = hitThres)
})

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
  
  updateAceEditor(session, editorId = "tableSumm_out", value = tableSummCode() )
  
  tableGroup
  
})

tableSummCode <- reactive({
  
  catType = as.numeric(input$radioMaxGroup)
  category <- c("Biological","Chemical","Chemical Class")[catType]
  hitThres <- hitThresValue()
  
  tableSummCode <- paste0(rCodeSetup(),"
# Use the table_tox_rank function for a formatted DT table
tableSum <- stats_of_groupings(chemicalSummary, 
                  category = '",category,"',
                  mean_logic = ",as.logical(input$meanEAR),",
                  hit_threshold = ",hitThres,")")
  
  return(tableSummCode)
  
})

output$downloadTable <- downloadHandler(
  filename = "tableSum.csv",
  content = function(file) {

    write.csv(tableSummData(), file = file, row.names = FALSE)
  }
)