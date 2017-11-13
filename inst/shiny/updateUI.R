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
                       group = init_Groups,
                       flags = initFlags,
                       sites = "All")

observeEvent(input$pickAssay, {
  epDF[["assays"]] <- NULL
  epDF[["assays"]] <- input$assay
  
})

observeEvent(input$pickFlags, {
  epDF[["flags"]] <- NULL
  epDF[["flags"]] <- input$flags
})

observeEvent(input$sites, {
  epDF[["sites"]] <- NULL
  epDF[["sites"]] <- input$sites
})

observeEvent(input$changeAnn, {
  epDF[["groupColName"]] <- NULL
  epDF[["groupColName"]] <- input$groupCol
  
})

observeEvent(input$allGroup, {
  epDF[["group"]] <- NULL
  epDF[["group"]] <- input$group
  
})

observeEvent(input$epGroup, {
  epDF[["epGroup"]] <- NULL
  epDF[["epGroup"]] <- input$epGroup
  
})

observe({
  
  groupCol <- epDF[["groupColName"]]
  assays <- epDF[["assays"]]
  
  ep <- filter_groups(ep = cleaned_ep, assays = assays, groupCol = groupCol, remove_groups = "")
  
  orderBy <- ep[,"groupCol"]
  orderNames <- names(table(orderBy))
  nEndPoints <- as.integer(table(orderBy))
  
  df <- data.frame(orderNames,nEndPoints,stringsAsFactors = FALSE) %>%
    arrange(desc(nEndPoints))
  
  dropDownHeader <- c(paste0(df$orderNames," (",df$nEndPoints,")"))
  
  selChoices <- df$orderNames
  
  if(epDF[["groupColName"]] == "intended_target_family"){
    selChoices <- selChoices[!(selChoices %in% c("Background Measurement","Undefined"))]
  }
  epDF[["group"]] <- NULL
  epDF[["group"]] <- selChoices
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

  chemicalSummary <- chemicalSummary()
  valueText <- "All"
  
  if (input$radioMaxGroup == 2){
    valueText <- c("All",unique(as.character(chemicalSummary$chnm)))
  } else if(input$radioMaxGroup == 3){
    valueText <- c("All",unique(as.character(chemicalSummary$Class)))
  } else if(input$radioMaxGroup == 1){
    valueText <- c("All",unique(as.character(chemicalSummary$Bio_category)))
  } 
  
  epDF[["epGroup"]] <- NULL
  
  if(length(valueText) > 1){
    epDF[["epGroup"]] <- valueText[2]
  } else {
    epDF[["epGroup"]] <- "All"
  }

  updateSelectInput(session, "epGroup", choices = valueText, selected = epDF[["epGroup"]])
})

