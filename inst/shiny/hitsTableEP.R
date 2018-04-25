output$hitsTableEPs <- DT::renderDataTable({
  
  validate(
    need(!is.null(input$data), "Please select a data set")
  )
  
  chemicalSummary <- chemicalSummary()
  catType <- as.numeric(input$radioMaxGroup)
  mean_logic <- input$meanEAR
  hitThres <- hitThresValue()
  
  tableEPs <- endpoint_hits_DT(chemicalSummary, 
                               category = c("Biological","Chemical","Chemical Class")[catType],
                               mean_logic = mean_logic,
                               hit_threshold = hitThres,
                               include_links = toxCast())
  updateAceEditor(session, editorId = "hitsTable_out", value = hitsTableEPCode() )
  tableEPs
})

hitsTableEPCode <- reactive({
  
  catType = as.numeric(input$radioMaxGroup)
  category <- c("Biological","Chemical","Chemical Class")[catType]
  hitThres <- hitThresValue()
  
  hitsTableEPCode <- paste0(rCodeSetup(),"
# Use the endpoint_hits_DT for a formatted DT table
hitTable <- endpoint_hits(chemicalSummary, 
              category = '",category,"',
              mean_logic = '",input$meanEAR,"',
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
  mean_logic <- input$meanEAR
  
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