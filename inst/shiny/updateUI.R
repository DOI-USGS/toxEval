choices <- reactive({
  rawData <- rawData()
  site_info <- rawData$chem_site
  choices <- c("All",site_info$`Short Name`)
  choices
})

observe({
  updateSelectInput(session, "sites", choices = choices())
})

epDF <- reactiveValues(assays = initAssay,
                       groupColName = "intended_target_family",
                       flags = flags)

observeEvent(input$pickAssay, {
  epDF[["assays"]] <- NULL
  epDF[["assays"]] <- input$assay
  
})

observeEvent(input$pickFlags, {
  epDF[["flags"]] <- NULL
  epDF[["flags"]] <- input$flags
})

observeEvent(input$changeAnn, {
  epDF[["groupColName"]] <- NULL
  epDF[["groupColName"]] <- input$groupCol
  
})



observe({
  
  groupCol <- epDF[["groupColName"]]
  assays <- epDF[["assays"]]
  
  ep <- data.frame(endPoint = cleaned_ep[["endPoint"]],
                   groupCol = cleaned_ep[[groupCol]],
                   assaysFull = cleaned_ep[["assay_source_name"]],
                   stringsAsFactors = FALSE) %>%
    filter(assaysFull %in% assays)
  
  orderBy <- ep[,"groupCol"]
  orderNames <- names(table(orderBy))
  nEndPoints <- as.integer(table(orderBy))
  
  df <- data.frame(orderNames,nEndPoints,stringsAsFactors = FALSE) %>%
    arrange(desc(nEndPoints))
  
  dropDownHeader <- c(paste0(df$orderNames," (",df$nEndPoints,")"))
  
  selChoices <- df$orderNames
  
  if(epDF[["groupColName"]] == "intended_target_family"){
    selChoices <- selChoices[!(selChoices %in% c("Background Measurement"))]
  }
  
  updateCheckboxGroupInput(session, "group", 
                           choices = setNames(df$orderNames,dropDownHeader),
                           selected = selChoices)
  
})

observe({
  labelText <- "Choose Chemical"
  
  if (input$radioMaxGroup == 2){
    labelText <- "Choose Chemical"
  } else if(input$radioMaxGroup == 3){
    labelText <- "Choose Class"
  } else if(input$radioMaxGroup == 1){
    labelText <- "Choose Group"
  }
  
  updateSelectInput(session, "epGroup", label = labelText)
})

observe({
  valueText <- "All"
  chemicalSummary <- chemicalSummary()
  
  if (input$radioMaxGroup == 2){
    uniqueChems <- c("All",unique(chemicalSummary$chnm))
    valueText <- uniqueChems
  } else if(input$radioMaxGroup == 3){
    valueText <- c("All",unique(chemicalSummary$Class))
  } else if(input$radioMaxGroup == 1){
    valueText <- c("All",unique(chemicalSummary$Bio_category))
  }
  
  updateSelectInput(session, "epGroup", choices = valueText, selected = valueText[2])
})