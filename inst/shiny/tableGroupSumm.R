output$tableGroupSumm <- DT::renderDataTable({
  
  validate(
    need(!is.null(input$data), "Please select a data set")
  )
  
  catType = as.numeric(input$radioMaxGroup)
  
  chemicalSummary <- chemicalSummary()
  hitThres <- hitThresValue()
  mean_logic <- as.logical(input$meanEAR)
  
  tableGroup <- table_tox_sum(chemicalSummary, 
                              category = c("Biological","Chemical","Chemical Class")[catType],
                              mean_logic = mean_logic,
                              hit_threshold = hitThres)
  
  updateAceEditor(session, editorId = "tableGroup_out", value = tableGroupCode() )
  tableGroup
  
})

output$nGroup <- renderUI({
  
  radio <- input$radioMaxGroup
  chemicalSummary <- chemicalSummary()
  siteToFind <- length(unique(chemicalSummary$site))
  
  if(radio == "2"){
    word <- "chemicals"
  } else if (radio == "3") {
    word <- "classes"
  } else if(radio == "4"){
    word <- "endPoints"
  } else {
    word <- "groups"
  }
  
  if(length(siteToFind) > 1){
    place <- "per site"
  } else {
    place <- ""
  }
  
  if(length(siteToFind) > 1){
    textUI <- paste("<h5>max = Maximum number of",word,"with hits per site</h5>",
                    "<h5>mean = Mean number of",word,"with hits per site</h5>",
                    "<h5>nSamples = Number of samples per site</h5>")
  } else {
    textUI <- paste("<h5>hits = Number of samples with hits </h5>",
                    "<h5>nSamples = Number of samples per site</h5>")
  }
  HTML(textUI)
  
})

tableGroupCode <- reactive({
  
  catType = as.numeric(input$radioMaxGroup)
  category <- c("Biological","Chemical","Chemical Class")[catType]
  hitThres <- hitThresValue()
  
  tableGroupCode <- paste0(rCodeSetup(),"
# Use the table_tox_sum function for a formatted DT table
tableGroupSum <- stats_of_hits(chemicalSummary, 
                  category = '",category,"',
                  mean_logic = ",as.logical(input$meanEAR),",
                  hit_threshold = ",hitThres,")")
  
  return(tableGroupCode)
  
})

tableSummGroupData <- reactive({
  
  validate(
    need(!is.null(input$data), "Please select a data set")
  )
  
  catType = as.numeric(input$radioMaxGroup)
  
  chemicalSummary <- chemicalSummary()
  hitThres <- hitThresValue()
  mean_logic <- as.logical(input$meanEAR)
  
  tableGroup <- stats_of_hits(chemicalSummary, 
                              category = c("Biological","Chemical","Chemical Class")[catType],
                              hit_threshold = hitThres)

})

output$downloadGroupTable <- downloadHandler(
  filename = "tableGroupSum.csv",
  content = function(file) {
    
    write.csv(tableSummGroupData(), file = file, row.names = FALSE)
  }
)