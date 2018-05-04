library(dplyr)
library(ggplot2)
library(DT)
library(leaflet)
library(data.table)
library(tidyr)
library(toxEval)
library(RColorBrewer)
library(grid)
library(RColorBrewer)
library(stringi)
library(readxl)

options(shiny.maxRequestSize=10*1024^2)

cleaned_ep <- clean_endPoint_info(endPointInfo) %>%
  mutate(endPoint = assay_component_endpoint_name)

choicesPerGroup <- apply(cleaned_ep, 2, function(x) length(unique(x[!is.na(x)])))
choicesPerGroup <- choicesPerGroup[which(as.numeric(choicesPerGroup) > 6)]
choicesPerGroup <- apply(cleaned_ep[,-1], 2, function(x) length(unique(x)))
groupChoices <- paste0(names(choicesPerGroup)," (",choicesPerGroup,")")

initAssay <- c("ATG","NVS","OT","TOX21","CEETOX", "APR", #"BSK"
               "CLD","TANGUAY","NHEERL_PADILLA",
               "NCCT_SIMMONS","ACEA")

init_Groups <- unique(cleaned_ep$intended_target_family)
init_Groups <- init_Groups[!is.na(init_Groups)]
init_Groups <- init_Groups[!(init_Groups %in% c("Background Measurement","Undefined"))]

all_flags <- c("Borderline",
                "OnlyHighest",
                "OneAbove",
                "Noisy",
                "HitCall",
                "GainAC50",
                "Biochemical")

initFlags <- c(#"Borderline",
                #"OnlyHighest",
                "OneAbove",
                "Noisy",
                "HitCall")
                #"GainAC50",
                #"Biochemical")

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

shinyServer(function(input, output,session) {
  
  observe({
    if (input$close > 0) stopApp()    
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
    setupCode <- paste0("######################################
# Setup:
library(toxEval)
#NOTE: Add path to path_to_file!!!
path_to_file <- '",fileName,"' 
tox_list <- create_toxEval(path_to_file)")
    
  if(toxCast()){ 
        setupCode <- paste0(setupCode,"
ACClong <- get_ACC(tox_list$chem_info$CAS)
ACClong <- remove_flags(ACClong = ACClong,
                        flagsShort = ",removeFlags,")

cleaned_ep <- clean_endPoint_info(endPointInfo)
filtered_ep <- filter_groups(cleaned_ep, 
                  groupCol = '",groupCol,"',
                  assays = ",assays,",
                  remove_groups = ",remove_groups,")

chemicalSummary <- get_chemical_summary(tox_list, ACClong, filtered_ep)")
  
  } else {
    setupCode <- paste0(setupCode,"
chemicalSummary <- get_chemical_summary(tox_list)" )  
  }
    
  if(sites != "All"){
    setupCode <- paste0(setupCode,"
site <- '",sites,"'
chemicalSummary <- chemicalSummary[chemicalSummary$shortName == site,]")
  }
    setupCode <- paste0(setupCode,"
######################################")
    return(setupCode)
  })
  
  chemicalSummary <- reactive({

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
      
      if(all(is.null(rawData$benchmarks))){

        ACClong <- get_ACC(rawData$chem_info$CAS)
        ACClong <- remove_flags(ACClong, flagsShort = removeFlags)
        
        remove_groups <- unique(cleaned_ep[[groupCol]])[which(!unique(cleaned_ep[[groupCol]]) %in% groups)]
        remove_groups <- remove_groups[!is.na(remove_groups)]
        
        filtered_ep <- filter_groups(cleaned_ep, 
                                     groupCol = groupCol, assays = assays,
                                     remove_groups = remove_groups)
        chemicalSummary <- get_chemical_summary(rawData,
                                                ACClong,
                                                filtered_ep) 
        toxCast_val <<- TRUE
        
      } else {
        chemicalSummary <- get_chemical_summary(rawData) 

        toxCast_val <<- FALSE
      }
      
    } else {
      chemicalSummary <- data.frame(casrn = character(),
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
    
    return(chemicalSummary)
    
  })
  
  toxCast <- reactive({
    rawData <- rawData()

    toxCast_val <- all(is.null(rawData$benchmarks))
    
    return(toxCast_val)
  })
  
  output$isTox <- reactive(toxCast())
  outputOptions(output, "isTox", suspendWhenHidden = FALSE)

  
  output$title_text <- renderText({
    
    if(toxCast()){
      textUI <- "Analysis using ToxCast endPoints"
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
                           "Chemical" = "for all chemicals",
                           "Biological" = "for chemicals within a grouping",
                           "Chemical Class" = "for chemicals within a class"
      )      
    }
    
    if(site == "All"){

      title <- paste("Summing EARs",pretty_cat, "for a given sample,")
      
      if(mean_logic){
        title <- paste(title,"taking the mean of each site")
      } else {
        title <- paste(title,"taking the max of each site")
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
      
      title <- paste(title,"
                     ", siteTable[["Fullname"]][which(siteTable$`Short Name` == site)])
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

###############################################################    
# Map Stuff:
  source("mapStuff.R",local=TRUE)$value
############################################################## 

###############################################################    
# Benchmark Stuff:
  source("benchmarks.R",local=TRUE)$value
############################################################## 


})