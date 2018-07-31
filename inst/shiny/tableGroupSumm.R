output$tableGroupSumm <- DT::renderDataTable({
  
  validate(
    need(!is.null(rawData_data$data), "Please select a data set")
  )
  
  catType = as.numeric(input$radioMaxGroup)
  
  chemical_summary <- chemical_summary()
  hitThres <- hitThresValue()
  sum_logic <- as.logical(input$sumEAR)
  
  tableGroup <- hits_summary_DT(chemical_summary, 
                              category = c("Biological","Chemical","Chemical Class")[catType],
                              sum_logic = sum_logic,
                              hit_threshold = hitThres)
  
  shinyAce::updateAceEditor(session, editorId = "tableGroup_out", value = tableGroupCode() )
  tableGroup
  
})

output$nGroup <- renderText({

  sum_logic <- as.logical(input$sumEAR)
  hit_thres <- hitThresValue()
  catType = as.numeric(input$radioMaxGroup)
  
  catType = as.numeric(input$radioMaxGroup)
  category <- c("group","chemical","chemical class")[catType]
  

  sum_word <- ifelse(sum_logic, "sum of EARs", "maximum EAR")

  textUI <- paste("Samples with hits = Number of samples where the",sum_word,"within a",category," is greater than",hit_thres)

  HTML(textUI)
  
})

tableGroupCode <- reactive({
  
  catType = as.numeric(input$radioMaxGroup)
  category <- c("Biological","Chemical","Chemical Class")[catType]
  hitThres <- hitThresValue()
  
  sum_logic <- as.logical(input$sumEAR)
  if(sum_logic){
    tableGroupCode <- paste0(rCodeSetup(),"
# Use the hits_summary_DT function for a formatted DT table
tableGroupSum <- hits_summary(chemical_summary, 
                  category = '",category,"',
                  hit_threshold = ",hitThres,")")    
  } else {
    tableGroupCode <- paste0(rCodeSetup(),"
# Use the hits_summary_DT function for a formatted DT table
                             tableGroupSum <- hits_summary(chemical_summary, 
                             category = '",category,"',
                             sum_logic = FALSE,
                             hit_threshold = ",hitThres,")")
  }

  
  return(tableGroupCode)
  
})

tableSummGroupData <- reactive({
  
  validate(
    need(!is.null(rawData_data$data), "Please select a data set")
  )
  
  catType = as.numeric(input$radioMaxGroup)
  
  chemical_summary <- chemical_summary()
  hitThres <- hitThresValue()
  sum_logic <- as.logical(input$sumEAR)
  
  tableGroup <- hits_summary(chemical_summary, 
                              category = c("Biological","Chemical","Chemical Class")[catType],
                              hit_threshold = hitThres,
                             sum_logic = sum_logic)

})

output$downloadGroupTable <- downloadHandler(
  filename = "tableGroupSum.csv",
  content = function(file) {
    
    write.csv(tableSummGroupData(), file = file, row.names = FALSE)
  }
)