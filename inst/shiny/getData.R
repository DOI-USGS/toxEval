rawData <- reactive({
  if(!is.null(input$data)){
    
    path <- file.path(input$data$datapath)
    newPath <- paste0(input$data$datapath,"_",input$data$name)
    newPath <- gsub(", ","_",newPath)
    file.rename(from = path, to = newPath)
    
    chem_data <- read_excel(newPath, sheet = "Data")
    chem_info <- read_excel(newPath, sheet = "Chemicals") 
    chem_site <- read_excel(newPath, sheet = "Sites")
    stationINFO <<- chem_site
    
    rawData <- list(chem_data=chem_data,
                    chem_info=chem_info,
                    chem_site=chem_site)
  } else {
    rawData <- NULL
  }
  
})

