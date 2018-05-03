output$hitsTable <- DT::renderDataTable({    
  validate(
    need(!is.null(input$data), "Please select a data set")
  )
  
  catType = as.numeric(input$radioMaxGroup)
  
  chemicalSummary <- chemicalSummary()
  hitThres <- hitThresValue()
  mean_logic <- as.logical(input$meanEAR)
  sum_logic <- as.logical(input$sumEAR)
  
  tableGroup <- hits_by_groupings_DT(chemicalSummary, 
                                   category = c("Biological","Chemical","Chemical Class")[catType],
                                   mean_logic = mean_logic,
                                   sum_logic = sum_logic,
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
  
  sum_logic <- as.logical(input$sumEAR)
  if(sum_logic){
    siteHitCode <- paste0(rCodeSetup(),"
# Use the hits_by_groupings_DT function for the formatted DT table
hitSiteTable <- hits_by_groupings(chemicalSummary, 
                    category = '",category,"',
                    mean_logic = ",input$meanEAR,",
                    hit_threshold = ",hitThres,")")
  } else {
    siteHitCode <- paste0(rCodeSetup(),"
# Use the hits_by_groupings_DT function for the formatted DT table
hitSiteTable <- hits_by_groupings(chemicalSummary, 
                    category = '",category,"',
                    mean_logic = ",input$meanEAR,",
                    sum_logic = FALSE,
                    hit_threshold = ",hitThres,")")    
  }
  return(siteHitCode)
  
})

siteHitTableData <- reactive({
  validate(
    need(!is.null(input$data), "Please select a data set")
  )
  
  catType = as.numeric(input$radioMaxGroup)
  
  chemicalSummary <- chemicalSummary()
  hitThres <- hitThresValue()
  mean_logic <- as.logical(input$meanEAR)
  sum_logic <- as.logical(input$sumEAR)
  
  tableGroup <- hits_by_groupings(chemicalSummary, 
                                   category = c("Biological","Chemical","Chemical Class")[catType],
                                   mean_logic = mean_logic,
                                  sum_logic = sum_logic,
                                   hit_threshold = hitThres)
})

output$downloadSiteHitTable <- downloadHandler(
  filename = "hitSiteTable.csv",
  content = function(file) {
    
    write.csv(siteHitTableData(), file = file, row.names = TRUE)
  }
)