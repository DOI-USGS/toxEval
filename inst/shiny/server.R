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
  rename(endPoint = assay_component_endpoint_name)

choicesPerGroup <- apply(cleaned_ep, 2, function(x) length(unique(x[!is.na(x)])))
choicesPerGroup <- which(cleaned_ep > 6 & cleaned_ep < 100)

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

choicesPerGroup <- apply(cleaned_ep[,-1], 2, function(x) length(unique(x)))
groupChoices <- paste0(names(choicesPerGroup)," (",choicesPerGroup,")")

initAssay <- c("ATG","NVS","OT","TOX21","CEETOX", "APR", #"BSK"
            "CLD","TANGUAY","NHEERL_PADILLA",
            "NCCT_SIMMONS","ACEA")

flags <- c( "OneAbove","Noisy","HitCall")


shinyServer(function(input, output,session) {
  
  observe({
    if (input$close > 0) stopApp()    
  })
  
  source("getData.R",local=TRUE)$value
  
  chemicalSummary <- reactive({
    
    groupCol <- epDF[["groupColName"]]
    assays <- epDF[["assays"]]
    flags <- epDF[["flags"]]
    sites <- epDF[["sites"]]
    
    rawData <- rawData()
    if(!is.null(rawData)){
      chem_data <- rawData$chem_data
      chem_info <- rawData$chem_info
      chem_site <- rawData$chem_site
      
      if(sites != "All"){
        chem_site <- chem_site[chem_site$`Short Name` == sites,]
        siteID <- chem_site$SiteID
        chem_data <- chem_data[chem_data$SiteID == siteID,]
        
      }
      
      ACClong <- get_ACC(chem_info$CAS)
      
      ACClong <- remove_flags(ACClong, flagsShort = epDF[["flags"]])
      
      remove_groups <- unique(cleaned_ep[[groupCol]])[which(!unique(cleaned_ep[[groupCol]]) %in% input$group)]
      # So maybe do this in clean?:
      cleaned_ep <- rename(cleaned_ep, assay_component_endpoint_name = endPoint)
      filtered_ep <- filter_groups(cleaned_ep, 
                                    groupCol = groupCol, 
                                    remove_groups = remove_groups)

      chemicalSummary <- get_chemical_summary(ACClong,
                                              filtered_ep,
                                              chem_data, 
                                              chem_site, 
                                              chem_info)  
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
    
  })
  
  hitThresValue <- eventReactive(input$changeHit, ignoreNULL = FALSE, {
    hitThresValue <- input$hitThres
    hitThresValue
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

  


})