rawData <- reactive({
  if(!is.null(input$data)){
    
    path <- file.path(input$data$datapath)
    
    newPath <- paste0(input$data$datapath,"_",input$data$name)
    
    epDF[["fileName"]] <- input$data$name
    
    newPath <- gsub(", ","_",newPath)
    
    file.rename(from = path, to = newPath)
    
    if(tools::file_ext(input$data$name) == "xlsx" |
       tools::file_ext(input$data$name) == "xls"){
      
      rawData <- create_toxEval(newPath)
       
      
    } else if(tools::file_ext(input$data$name) == "RData" | 
              tools::file_ext(input$data$name) == "rda"){
      load(newPath)
      
      if(exists("chem_data") & exists("chem_info") & 
         exists("chem_site")) {
        if(exists("exlusions")){
                  rawData <- create_toxEval(chem_data=chem_data,
                        chem_info=chem_info,
                        chem_site=chem_site,
                        exclusions=exclusions)
        } else {
          rawData <- create_toxEval(chem_data=chem_data,
                          chem_info=chem_info,
                          chem_site=chem_site,
                          exclusions=NULL)
        }
      } 

      if(any(sapply(ls(), function(x) class(get(x))) == "toxEval")){
        rawData <- get(ls()[which(sapply(ls(), function(x) class(get(x))) == "toxEval")])
      }
      
    } else if(tools::file_ext(input$data$name) == "rds"){
      rawData <- readRDS(newPath)
    }


  } else {
    rawData <- NULL
  }
  
})

