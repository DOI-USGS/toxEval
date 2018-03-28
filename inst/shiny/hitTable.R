output$hitsTable <- DT::renderDataTable({    
  validate(
    need(!is.null(input$data), "Please select a data set")
  )
  
  catType = as.numeric(input$radioMaxGroup)
  
  chemicalSummary <- chemicalSummary()
  hitThres <- hitThresValue()
  mean_logic <- as.logical(input$meanEAR)
  
  tableGroup <- table_tox_endpoint(chemicalSummary, 
                                   category = c("Biological","Chemical","Chemical Class")[catType],
                                   mean_logic = mean_logic,
                                   hit_threshold = hitThres)
  
  updateAceEditor(session, editorId = "siteHit_out", value = siteHitCode() )
  tableGroup
  
})

output$siteHitText <- renderUI({
  
  if(input$sites == "All"){
    HTML(paste("<h4>Number of sites with hits</h4>"))
  } else {
    HTML(paste("<h4>Number of samples with hits</h4>"))
  }
  
})

siteHitCode <- reactive({
  
  catType = as.numeric(input$radioMaxGroup)
  category <- c("Biological","Chemical","Chemical Class")[catType]
  hitThres <- hitThresValue()
  
  siteHitCode <- paste0(rCodeSetup(),"
table_tox_endpoint(chemicalSummary, 
                    category = '",category,"',
                    mean_logic = ",as.logical(input$meanEAR),",
                    hit_threshold = ",hitThres,")")
  
  return(siteHitCode)
  
})