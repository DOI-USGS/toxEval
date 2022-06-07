tableSummData <- reactive({
  catType = as.numeric(input$radioMaxGroup)
  
  chemical_summary <- chemical_summary()
  hitThres <- hitThresValue()
  mean_logic <- as.logical(input$meanEAR)
  sum_logic <- as.logical(input$sumEAR)

  tableGroup <- rank_sites(chemical_summary, 
                           category = c("Biological","Chemical","Chemical Class")[catType],
                           mean_logic = mean_logic,
                           sum_logic = sum_logic,
                           hit_threshold = hitThres)
})

output$tableSumm <- DT::renderDataTable({
  
  validate(
    need(!is.null(rawData_data$data), "")
  )
  
  catType = as.numeric(input$radioMaxGroup)
  
  chemical_summary <- chemical_summary()
  hitThres <- hitThresValue()
  mean_logic <- as.logical(input$meanEAR)
  sum_logic <- as.logical(input$sumEAR)
  
  tableGroup <- rank_sites_DT(chemical_summary, 
                               category = c("Biological","Chemical","Chemical Class")[catType],
                               mean_logic = mean_logic,
                              sum_logic = sum_logic,
                               hit_threshold = hitThres)
  
  shinyAce::updateAceEditor(session, editorId = "tableSumm_out", value = tableSummCode() )
  
  tableGroup
  
}, server = FALSE)

tableSummCode <- reactive({
  
  catType = as.numeric(input$radioMaxGroup)
  category <- c("Biological","Chemical","Chemical Class")[catType]
  hitThres <- hitThresValue()
  sum_logic <- as.logical(input$sumEAR)
  
  if(sum_logic){
  tableSummCode <- paste0(rCodeSetup(),"
# Use the rank_sites_DT function for a formatted DT table
tableSum <- rank_sites(chemical_summary, 
                  category = '",category,"',
                  mean_logic = ",as.logical(input$meanEAR),",
                  hit_threshold = ",hitThres,")")
  } else {
    tableSummCode <- paste0(rCodeSetup(),"
# Use the rank_sites_DT function for a formatted DT table
tableSum <- rank_sites(chemical_summary, 
                      category = '",category,"',
                      mean_logic = ",as.logical(input$meanEAR),",
                      sum_logic = FALSE,
                      hit_threshold = ",hitThres,")")    
  }
  return(tableSummCode)
  
})

output$downloadTable <- downloadHandler(
  filename = "tableSum.csv",
  content = function(file) {

    write.csv(tableSummData(), file = file, row.names = FALSE)
  }
)