tableSummData <- reactive({
  catType = as.numeric(input$radioMaxGroup)
  
  chemicalSummary <- chemicalSummary()
  hitThres <- hitThresValue()
  mean_logic <- input$meanEAR
  
  tableGroup <- rank_sites(chemicalSummary, 
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
  mean_logic <- input$meanEAR
  
  tableGroup <- rank_sites_DT(chemicalSummary, 
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
# Use the rank_sites_DT function for a formatted DT table
tableSum <- rank_sites(chemicalSummary, 
                  category = '",category,"',
                  mean_logic = '",input$meanEAR,"',
                  hit_threshold = ",hitThres,")")
  
  return(tableSummCode)
  
})

output$downloadTable <- downloadHandler(
  filename = "tableSum.csv",
  content = function(file) {

    write.csv(tableSummData(), file = file, row.names = FALSE)
  }
)