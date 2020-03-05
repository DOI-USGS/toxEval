output$hitsTableEPs <- DT::renderDataTable({
  
  validate(
    need(!is.null(rawData_data$data), "")
  )
  
  chemical_summary <- chemical_summary()
  catType <- as.numeric(input$radioMaxGroup)
  mean_logic <- as.logical(input$meanEAR)
  sum_logic <- as.logical(input$sumEAR)
  hitThres <- hitThresValue()

  tableEPs <- endpoint_hits_DT(chemical_summary, 
                               category = c("Biological","Chemical","Chemical Class")[catType],
                               mean_logic = mean_logic,
                               sum_logic = sum_logic,
                               hit_threshold = hitThres,
                               include_links = toxCast())
  shinyAce::updateAceEditor(session, editorId = "hitsTable_out", value = hitsTableEPCode() )
  tableEPs
})

hitsTableEPCode <- reactive({
  
  catType = as.numeric(input$radioMaxGroup)
  category <- c("Biological","Chemical","Chemical Class")[catType]
  hitThres <- hitThresValue()
  sum_logic <- input$sumEAR
  if(sum_logic){
    hitsTableEPCode <- paste0(rCodeSetup(),"
# Use the endpoint_hits_DT for a formatted DT table
hitTable <- endpoint_hits(chemical_summary, 
              category = '",category,"',
              mean_logic = ",as.logical(input$meanEAR),",
              hit_threshold = ",hitThres,")")
  } else {
    hitsTableEPCode <- paste0(rCodeSetup(),"
# Use the endpoint_hits_DT for a formatted DT table
hitTable <- endpoint_hits(chemical_summary, 
              category = '",category,"',
              mean_logic = ",as.logical(input$meanEAR),",
              sum_logic = FALSE,
              hit_threshold = ",hitThres,")")    
  }
  return(hitsTableEPCode)
  
})

hitTableData <- reactive({
  validate(
    need(!is.null(rawData_data$data), "Please select a data set")
  )
  
  catType = as.numeric(input$radioMaxGroup)
  
  chemical_summary <- chemical_summary()
  hitThres <- hitThresValue()
  mean_logic <- as.logical(input$meanEAR)
  sum_logic <- as.logical(input$sumEAR)
  tableGroup <- endpoint_hits(chemical_summary, 
                             category = c("Biological","Chemical","Chemical Class")[catType],
                             mean_logic = mean_logic,
                             sum_logic = sum_logic,
                             hit_threshold = hitThres)
})

output$downloadHitTable <- downloadHandler(
  filename = "hitTable.csv",
  content = function(file) {
    
    write.csv(hitTableData(), file = file, row.names = FALSE)
  }
)

output$epHitTitle <- renderText({
  
  validate(
    need(!is.null(rawData_data$data), "Please select a data set")
  )
  
  hitThres <- hitThresValue()
  mean_logic <- as.logical(input$meanEAR)
  sum_logic <- as.logical(input$sumEAR)  
  
  site_word <- ifelse(input$sites == "All","sites","samples")
  
  sum_word <- ifelse(sum_logic, "sum","max")
  
  catType = as.numeric(input$radioMaxGroup)
  category <- c("group","chemical","chemical class")[catType]
  
  text_ui <- paste("Number of",site_word,"where the",sum_word,
                   "of the EARs in a", category,"is greater than",hitThres)
  
  
  return(HTML(text_ui))
  
})