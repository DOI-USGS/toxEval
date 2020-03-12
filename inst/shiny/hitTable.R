output$hitsTable <- DT::renderDataTable({    
  validate(
    need(!is.null(rawData_data$data), "")
  )
  
  catType = as.numeric(input$radioMaxGroup)
  
  chemical_summary <- chemical_summary()
  hitThres <- hitThresValue()
  mean_logic <- as.logical(input$meanEAR)
  sum_logic <- as.logical(input$sumEAR)
  
  tableGroup <- hits_by_groupings_DT(chemical_summary, 
                                   category = c("Biological","Chemical","Chemical Class")[catType],
                                   mean_logic = mean_logic,
                                   sum_logic = sum_logic,
                                   hit_threshold = hitThres)
  
  shinyAce::updateAceEditor(session, editorId = "siteHit_out", value = siteHitCode() )
  tableGroup
  
})

output$siteHitText <- renderText({
  
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

siteHitCode <- reactive({
  
  catType = as.numeric(input$radioMaxGroup)
  category <- c("Biological","Chemical","Chemical Class")[catType]
  hitThres <- hitThresValue()
  
  sum_logic <- as.logical(input$sumEAR)
  if(sum_logic){
    siteHitCode <- paste0(rCodeSetup(),"
# Use the hits_by_groupings_DT function for the formatted DT table
hitSiteTable <- hits_by_groupings(chemical_summary, 
                    category = '",category,"',
                    mean_logic = ",input$meanEAR,",
                    hit_threshold = ",hitThres,")")
  } else {
    siteHitCode <- paste0(rCodeSetup(),"
# Use the hits_by_groupings_DT function for the formatted DT table
hitSiteTable <- hits_by_groupings(chemical_summary, 
                    category = '",category,"',
                    mean_logic = ",input$meanEAR,",
                    sum_logic = FALSE,
                    hit_threshold = ",hitThres,")")    
  }
  return(siteHitCode)
  
})

siteHitTableData <- reactive({
  validate(
    need(!is.null(rawData_data$data), "Please select a data set")
  )
  
  catType = as.numeric(input$radioMaxGroup)
  
  chemical_summary <- chemical_summary()
  hitThres <- hitThresValue()
  mean_logic <- as.logical(input$meanEAR)
  sum_logic <- as.logical(input$sumEAR)
  
  tableGroup <- hits_by_groupings(chemical_summary, 
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