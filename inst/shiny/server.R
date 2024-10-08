
library(leaflet)
library(shiny)

options(shiny.maxRequestSize=70*1024^2)

cleaned_ep <- clean_endPoint_info(end_point_info) %>%
  dplyr::mutate(endPoint = assay_component_endpoint_name)

choicesPerGroup <- apply(cleaned_ep, 2, function(x) length(unique(x[!is.na(x)])))
choicesPerGroup <- choicesPerGroup[which(as.numeric(choicesPerGroup) > 6)]
choicesPerGroup <- apply(cleaned_ep[,-1], 2, function(x) length(unique(x)))
groupChoices <- paste0(names(choicesPerGroup)," (",choicesPerGroup,")")

initAssay <- c("BSK")

init_Groups <- unique(cleaned_ep$intended_target_family)
init_Groups <- init_Groups[!is.na(init_Groups)]
init_Groups <- init_Groups[!(init_Groups %in% c("Background Measurement","Undefined"))]

all_flags <- flags$flag_id

initFlags <-  c(5, 6, 11, 15, 18)

sitesOrdered <- c("StLouis","Pigeon","Nemadji","WhiteWI","Bad","Montreal","PresqueIsle",
                  "Ontonagon","Sturgeon","Tahquamenon",
                  "Burns","IndianaHC","StJoseph","PawPaw",        
                  "Kalamazoo","GrandMI","Milwaukee","Muskegon",      
                  "WhiteMI","Sheboygan","PereMarquette","Manitowoc",    
                  "Manistee","Fox","Oconto","Peshtigo",      
                  "Menominee","Indian","Cheboygan","Ford",         
                  "Escanaba","Manistique",
                  "ThunderBay","AuSable","Rifle",
                  "Saginaw","BlackMI","Clinton","Rouge","HuronMI","Raisin","Maumee",
                  "Portage","Sandusky","HuronOH","Vermilion","BlackOH","Rocky","Cuyahoga","GrandOH",
                  "Cattaraugus","Tonawanda","Genesee","Oswego","BlackNY","Oswegatchie","Grass","Raquette","StRegis")

great_lakes <- c("Lake Superior",
                 "Lake Michigan",
                 "Lake Huron",
                 "Lake Erie",
                 "Lake Ontario")

toxCast_val <<- TRUE
latest_map <<- NULL

shiny::shinyServer(function(input, output,session) {
  
  observe({
    if (input$close > 0) shiny::stopApp()    
  })

  source("getData.R",local=TRUE)$value
  
  rCodeSetup <- reactive({
    groupCol <- epDF[["groupColName"]]
    assays <- epDF[["assays"]]
    flags <- epDF[["flags"]]
    sites <- epDF[["sites"]]
    groups <- epDF[["group"]]
    fileName <- epDF[["fileName"]]

    sites <- input$sites
    
    remove_groups <- unique(cleaned_ep[[groupCol]])[which(!unique(cleaned_ep[[groupCol]]) %in% groups)]
    remove_groups <- remove_groups[!is.na(remove_groups)]
    
    removeFlags <- all_flags[!(all_flags %in% flags)]
    
    removeFlags <- paste0("c('",paste0(removeFlags, collapse = "','"),"')")
    assays <- paste0("c('",paste0(assays, collapse = "','"),"')")
    remove_groups <- paste0("c('",paste0(remove_groups, collapse = "','"),"')")
    setupCode <- paste0("#### Setup ####
library(toxEval)
#NOTE: Add path to path_to_file!!!
path_to_file <- '",fileName,"' 
tox_list <- create_toxEval(path_to_file)")
    
  if(toxCast()){ 
        setupCode <- paste0(setupCode,"
ACC <- get_ACC(tox_list$chem_info$CAS)
ACC <- remove_flags(ACC = ACC,
                    flag_id = ",removeFlags,")

cleaned_ep <- clean_endPoint_info(end_point_info)
filtered_ep <- filter_groups(cleaned_ep, 
                  groupCol = '",groupCol,"',
                  remove_assays = ",assays,",
                  remove_groups = ",remove_groups,")

chemical_summary <- get_chemical_summary(tox_list, ACC, filtered_ep)")
  
  } else {
    setupCode <- paste0(setupCode,"
chemical_summary <- get_chemical_summary(tox_list)" )  
  }
    
  if(sites != "All"){
    setupCode <- paste0(setupCode,"
site <- '",sites,"'
chemical_summary <- chemical_summary[chemical_summary$shortName == site,]")
  }
    setupCode <- paste0(setupCode,"
######################################")
    return(setupCode)
  })
  
  chemical_summary <- reactive({

    validate(
      need(!is.null(rawData_data$data), "Please select a data set")
    )
    
    groupCol <- epDF[["groupColName"]]
    assays <- epDF[["assays"]]
    flags <- epDF[["flags"]]
    sites <- epDF[["sites"]]
    groups <- epDF[["group"]]
    
    removeFlags <- all_flags[!(all_flags %in% flags)]

    rawData <- rawData()

    if(!is.null(rawData)){

      if(sites != "All"){
        rawData$chem_site <- rawData$chem_site[rawData$chem_site$`Short Name` == sites,]
        siteID <- rawData$chem_site$SiteID
        rawData$chem_data <- rawData$chem_data[rawData$chem_data$SiteID == siteID,]
        
      }

      if(all(is.null(rawData$benchmarks)) || 
         nrow(rawData$benchmarks) == 0 ||
         as.logical(input$useToxCast)){

        ACC <- get_ACC(rawData$chem_info$CAS)
        ACC <- remove_flags(ACC, flag_id = removeFlags)
        
        remove_groups <- unique(cleaned_ep[[groupCol]])[which(!unique(cleaned_ep[[groupCol]]) %in% groups)]
        remove_groups <- remove_groups[!is.na(remove_groups)]
        
        filtered_ep <- filter_groups(cleaned_ep, 
                                     groupCol = groupCol, 
                                     remove_assays = assays,
                                     remove_groups = remove_groups)
        chemical_summary <- get_chemical_summary(rawData,
                                                ACC,
                                                filtered_ep) 
        toxCast_val <<- TRUE
        
      } else {
        chemical_summary <- get_chemical_summary(rawData) 

        toxCast_val <<- FALSE
      }
      
    } else {
      chemical_summary <- data.frame(casrn = character(),
                       chnm = character(),
                       endPoint = character(),
                       site = character(),
                       date = numeric(),
                       EAR = numeric(),
                       shortName = character(),
                       Class = character(),
                       Bio_category = character(),
                       stringsAsFactors = FALSE)
    }
    validate(
      need(nrow(chemical_summary) > 0, "No data matched this combination")
    )
    
    return(chemical_summary)
    
  })
  
  toxCast <- reactive({
    rawData <- rawData()
    
    if(all(is.null(rawData$benchmarks))){
      toxCast_val <- TRUE
    } else {
      toxCast_val <- as.logical(input$useToxCast)
    }
    
    
    return(toxCast_val)
  })
  
  hasBenchmarks <- reactive({
    rawData <- rawData()
    return(!all(is.null(rawData$benchmarks)))
  })
  
  output$isTox <- reactive(toxCast())
  output$hasBenchmarks <- reactive(hasBenchmarks())
  outputOptions(output, "isTox", suspendWhenHidden = FALSE)
  outputOptions(output, "hasBenchmarks", suspendWhenHidden = FALSE)

  output$title_text <- renderText({
    
    if(toxCast()){
      textUI <- paste("Analysis using ToxCast", toxEval:::dbVersion())
    } else {
      textUI <- "Analysis using CUSTOM endPoints:
      Many dropdowns on sidebar will have no effect"
    }
    
    HTML(textUI)
    })
  
  output$siteText <- renderText({

    site <- input$sites
    
    if(site == "All"){
      textUI <- ""
    } else {
      siteTable <- rawData()[["chem_site"]]
      textUI <- siteTable[["Fullname"]][which(siteTable$`Short Name` == site)]
    }

    HTML(textUI)

  })
  
  output$meanText <- renderText({
    validate(
      need(!is.null(rawData_data$data), "")
    )
    
    catType = as.numeric(input$radioMaxGroup)
    category <- c("group","chemical","chemical class")[catType]
    
    mean_logic <- as.logical(input$meanEAR)
    sum_logic <- as.logical(input$sumEAR)
    
    mean_word <- ifelse(mean_logic,"mean","max")
    sum_word <- ifelse(sum_logic,"Summation of EARs","Maximum EAR")
    textUI <- paste0(mean_word,"EAR = ",sum_word," within a ", category," per sample, ",mean_word," at each site")
    
    HTML(textUI)
  })
  
  output$freqText <- renderText({
    
    validate(
      need(!is.null(rawData_data$data), "Please select a data set")
    )
    
    catType = as.numeric(input$radioMaxGroup)
    category <- c("group","chemical","chemical class")[catType]
    
    hit_thres <- hitThresValue()
    sum_logic <- as.logical(input$sumEAR)
    if(sum_logic){
      textUI <- paste("freq = Fraction of samples where the sum of EARs within a specified",category,"is greater than",hit_thres)
    } else {
      textUI <- paste("freq = Fraction of samples where the max EAR within a specified",category,"is greater than",hit_thres)      
    }
    
    HTML(textUI)
  })
  
  output$columnText <- renderText({
    
    if(toxCast()){
      textUI <- "Analysis using ToxCast endPoints"
    } else {
      textUI <- "Analysis using CUSTOM endPoints:
      Many dropdowns on sidebar will have no effect"
    }
    
    HTML(textUI)
  })
  
  hitThresValue <- eventReactive(input$changeHit, ignoreNULL = FALSE, {
    hitThresValue <- input$hitThres
    hitThresValue
  })
  
  genericTitle <- reactive({

    tab <- input$mainOut
    
    catType = as.numeric(input$radioMaxGroup)
    category <- c("Biological","Chemical","Chemical Class")[catType]
    
    site <- input$sites
    siteTable <- rawData()[["chem_site"]]
    
    mean_logic <- as.logical(input$meanEAR)
    sum_logic <- as.logical(input$sumEAR)

    if(tab == "endpoint"){
      filterBy <- epDF[['epGroup']]
      
      pretty_cat <- switch(category, 
                           "Chemical" = paste("for",filterBy),
                           "Biological" = paste("for chemicals within the",filterBy,"class"),
                           "Chemical Class" = paste("for chemicals within the",filterBy,"group")
      )
    } else {
      pretty_cat <- switch(category, 
                           "Chemical" = "",
                           "Biological" = "for chemicals within a grouping",
                           "Chemical Class" = "for chemicals within a class"
      )      
    }
    
    if(site == "All"){
      
      if(sum_logic){
        title <- paste("Summing EARs",pretty_cat, "of a sample,")
      } else {
        title <- paste("Max EARs",pretty_cat, "of a sample,")
      }
      
      if(mean_logic){
        title <- paste(title,"\ntaking the mean of each site")
      } else {
        title <- paste(title,"\ntaking the max of each site")
      }
    } else {
      
      if(tab == "endpoint"){
        filterBy <- epDF[['epGroup']]
        pretty_cat <- switch(category, 
                             "Chemical" = filterBy,
                             "Biological" = paste("chemicals within",filterBy),
                             "Chemical Class" = paste("chemicals within",filterBy)
        )
      } else {
        pretty_cat <- switch(category, 
                             "Chemical" = "chemical",
                             "Biological" = "grouping",
                             "Chemical Class" = "chemical class"
        )
      }
      
      title <- paste("EAR per",pretty_cat)
      
      if(tab == "summaryBar"){
        title <- paste(title, "for individual samples")
      }
      
      title <- paste(title,"\n", siteTable[["Fullname"]][which(siteTable$`Short Name` == site)])
    }
    return(title)
    
  })
  
  source("updateUI.R",local=TRUE)$value
  
#############################################################   
# DT tables:
  source("tableGroupSumm.R",local=TRUE)$value
  source("tableSum.R",local=TRUE)$value
  source("hitTable.R",local=TRUE)$value
  source("hitsTableEP.R",local=TRUE)$value
###################################################################
 
###############################################################
# Graphs:   
  source("stackPlot.R",local=TRUE)$value
  source("boxPlot.R",local=TRUE)$value
  source("heatMap.R",local=TRUE)$value
  source("endpointGraph.R",local=TRUE)$value
################################################################
  output$mymap <- renderLeaflet({
    leaflet() %>%
      addProviderTiles("CartoDB.Positron") %>%
      setView(lng = -83.5, lat = 44.5, zoom=6)
  })
  
###############################################################    
  source("mapStuff.R",local=TRUE)$value
############################################################## 

###############################################################    
# Benchmark Stuff:
  source("benchmarks.R",local=TRUE)$value
############################################################## 


  session$onSessionEnded(stopApp)
})