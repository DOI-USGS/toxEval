rawData <- reactive({
  if(!is.null(input$data)){
    
    path <- file.path(input$data$datapath)
    
    newPath <- paste0(input$data$datapath,"_",input$data$name)
    newPath <- gsub(", ","_",newPath)
    
    file.rename(from = path, to = newPath)
    
    if(tools::file_ext(input$data$name) == "xlsx" |
       tools::file_ext(input$data$name) == "xls"){
      
      rawData <- create_toxEval(newPath)
       
      
    } else if(tools::file_ext(input$data$name) == "RData" | 
              tools::file_ext(input$data$name) == "rds"){
      load(newPath)
      
      rawData <- list(chem_data=chem_data,
                      chem_info=chem_info,
                      chem_site=chem_site,
                      exclusions=exclusions)
    }


  } else {
    rawData <- NULL
  }
  
})

