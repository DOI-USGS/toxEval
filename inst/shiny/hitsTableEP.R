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
  updateAceEditor(session, editorId = "hitsTable_out", value = hitsTableEPCode() )
  tableEPs
})

hitsTableEPCode <- reactive({
  
  catType = as.numeric(input$radioMaxGroup)
  category <- c("Biological","Chemical","Chemical Class")[catType]
  hitThres <- hitThresValue()
  
  hitsTableEPCode <- paste0(rCodeSetup(),"
# Use the table_endpoint_hits for a formatted DT table
hitTable <- endpointHits(chemicalSummary, 
              category = '",category,"',
              mean_logic = ",as.logical(input$meanEAR),",
              hit_threshold = ",hitThres,")")
  
  return(hitsTableEPCode)
  
})

hitTableData <- reactive({
  validate(
    need(!is.null(input$data), "Please select a data set")
  )
  
  catType = as.numeric(input$radioMaxGroup)
  
  chemicalSummary <- chemicalSummary()
  hitThres <- hitThresValue()
  mean_logic <- as.logical(input$meanEAR)
  
  tableGroup <- endpointHits(chemicalSummary, 
                             category = c("Biological","Chemical","Chemical Class")[catType],
                             mean_logic = mean_logic,
                             hit_threshold = hitThres)
})

output$downloadHitTable <- downloadHandler(
  filename = "hitTable.csv",
  content = function(file) {
    
    write.csv(hitTableData(), file = file, row.names = FALSE)
  }
)