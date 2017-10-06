rawData <- reactive({
  if(!is.null(input$data)){
    
    path <- file.path(input$data$datapath)
    
    newPath <- paste0(input$data$datapath,"_",input$data$name)
    newPath <- gsub(", ","_",newPath)
    
    file.rename(from = path, to = newPath)
    
    if(tools::file_ext(input$data$name) == "xlsx" |
       tools::file_ext(input$data$name) == "xls"){
      chem_data <- read_excel(newPath, sheet = "Data")
      chem_info <- read_excel(newPath, sheet = "Chemicals") 
      chem_site <- read_excel(newPath, sheet = "Sites") 
      if("Exclude" %in% excel_sheets(newPath)){
        exclusions <- read_excel(newPath, sheet = "Exclude") 
      } else {
        exclusions <- NULL
      }
       
      
    } else if(tools::file_ext(input$data$name) == "RData" | 
              tools::file_ext(input$data$name) == "rds"){
      load(newPath)
    }

    rawData <- list(chem_data=chem_data,
                    chem_info=chem_info,
                    chem_site=chem_site,
                    exclusions=exclusions)
  } else {
    rawData <- NULL
  }
  
})

