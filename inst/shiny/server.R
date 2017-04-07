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

createLink <- function(cas,ep, hits) {
  paste0('<a href="http://actor.epa.gov/dashboard/#chemical/',cas,'" target="_blank" >',hits,'</a>')
}

shinyServer(function(input, output,session) {
  
  source("getData.R",local=TRUE)$value
  
  chemicalSummary <- reactive({
    
    groupCol <- epDF[["groupColName"]]
    assays <- epDF[["assays"]]
    flags <- epDF[["flags"]]
    
    rawData <- rawData()
    if(!is.null(rawData)){
      chem_data <- rawData$chem_data
      chem_info <- rawData$chem_info
      chem_site <- rawData$chem_site
      
      #Remove DEET:
      chem_data <- chem_data[chem_data$CAS != "134-62-3",] 
      chem_info <- chem_info[chem_info$CAS != "134-62-3",] 
      
      ACClong <- get_ACC(chem_info$CAS)
      
      ACClong <- remove_flags(ACClong, flagsShort = epDF[["flags"]])
      
      remove_groups <- unique(cleaned_ep[[groupCol]])[which(!unique(cleaned_ep[[groupCol]]) %in% input$group)]
      # So maybe do this in clean?:
      cleaned_ep <- rename(cleaned_ep, assay_component_endpoint_name = endPoint)
      filtered_ep <<- filter_groups(cleaned_ep, 
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
  
  hitThresValue <- eventReactive(input$changeHit, ignoreNULL = FALSE, {
    hitThresValue <- input$hitThres
    hitThresValue
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
  
#############################################################   
# DT tables:
  source("tableGroupSumm.R",local=TRUE)$value
  source("tableSum.R",local=TRUE)$value
  source("hitTable.R",local=TRUE)$value
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


#   
#   output$hitsTableEPs <- DT::renderDataTable({
# 
#     boxData <- chemicalSummary()
#     meanEARlogic <- as.logical(input$meanEAR)
#     catType <- as.numeric(input$radioMaxGroup)
#     
#     hitThres <- hitThresValue()
# 
#     fullData_init <- data.frame(Endpoint="",stringsAsFactors = FALSE)
#     fullData <- fullData_init
#     
#     if(length(unique(boxData$site)) > 1){
#       
#       for(i in unique(boxData$choices)){
#         dataSub <- boxData %>%
#           filter(choices == i) %>%
#           group_by(site, category, endPoint, date) %>%
#           summarize(sumEAR = sum(EAR)) %>%
#           group_by(site, category, endPoint) %>%
#           summarize(meanEAR = ifelse(meanEARlogic, mean(sumEAR),max(sumEAR))) %>%
#           group_by(category, endPoint) %>%
#           summarize(nSites = sum(meanEAR>hitThres)) %>%
#           data.frame() %>%
#           arrange(desc(nSites)) %>%
#           reshape(idvar="endPoint",timevar="category", direction="wide") 
#         
#         names(dataSub) <- gsub("nSites.","",names(dataSub))
#         names(dataSub)[1] <- "Endpoint"
# 
#         if(ncol(dataSub) > 2){
#           dataSub <- dataSub[,c(1,1+which(colSums(dataSub[,-1],na.rm = TRUE) != 0))]
#         }
#         
#         if(is.data.frame(dataSub)){
#           if(ncol(dataSub) > 2){
#             dataSub <- dataSub[(rowSums(dataSub[,-1],na.rm = TRUE) != 0),]
#           } else {
#             dataSub <- dataSub[which(dataSub[,-1] != 0 ),]
#           }
#           
#           dataSub <- dataSub %>%
#             data.frame() %>%
#             mutate(Group = i)
#           
#           fullData <- full_join(fullData,dataSub)
# 
#         }
#       }
#       
#     } else {
#       
#       for(i in unique(boxData$choices)){
#         dataSub <- boxData %>%
#           filter(choices == i) %>%
#           group_by(category, endPoint, date) %>%
#           summarise(sumEAR=sum(EAR)) %>%
#           data.frame() %>%
#           group_by(endPoint, category) %>%
#           summarise(nSites = sum(sumEAR>hitThres))%>%
#           data.frame() %>%
#           arrange(desc(nSites)) %>%
#           reshape(idvar="endPoint",timevar="category", direction="wide") 
#         
#         names(dataSub) <- gsub("nSites.","",names(dataSub))
#         names(dataSub)[1] <- "Endpoint"
#         
#         if(ncol(dataSub) > 2){
#           dataSub <- dataSub[,c(1,1+which(colSums(dataSub[,-1],na.rm = TRUE) != 0))]
#         }
#         
#         if(is.data.frame(dataSub)){
#           if(ncol(dataSub) > 2){
#             dataSub <- dataSub[(rowSums(dataSub[,-1],na.rm = TRUE) != 0),]
#           } else {
#             dataSub <- dataSub[which(dataSub[,-1] != 0 ),]
#           }
#           dataSub <- dataSub %>%
#             data.frame() %>%
#             mutate(Group = i)
#           
#           fullData <- full_join(fullData,dataSub)
#         }
#       }
#     }
#     
#     fullData <- fullData[,c("Endpoint","Group",names(fullData)[!(names(fullData) %in% c("Endpoint","Group"))])]
# 
#     fullData <- fullData[-1,]
#     
#     names(fullData) <- gsub("\\."," ",names(fullData))
#     names(fullData)[names(fullData) == "Human Drug  Non Prescription"] <- "Human Drug"
#     names(fullData)[names(fullData) == "Flavor/Fragrance"] <- "Flavor / Fragrance"
#     
#     sumOfColumns <- colSums(fullData[c(-1,-2)],na.rm = TRUE)
#     orderData <- order(sumOfColumns,decreasing = TRUE) 
#     orderData <- orderData[sumOfColumns[orderData] != 0] + 2
#     
#     fullData <- fullData[,c(1,2,orderData)]
#     colors <- brewer.pal(9,"Blues") #"RdYlBu"
#     
#     groups <- fullData$Groups
#     
#     if(catType == 2){
#       casKey <- select(boxData, chnm, casrn) %>%
#         distinct()
#       
#       hits <- sapply(fullData, function(x) as.character(x))
#       
#       for(k in 1:nrow(fullData)){
#         for(z in 3:ncol(fullData)){
#           if(!is.na(fullData[k,z])){
#             hits[k,z] <- createLink(ep = fullData$Endpoint[k], 
#                                     cas = casKey$casrn[casKey$chnm == names(fullData)[z]], 
#                                     hits = fullData[k,z])
#           }
#         }
#       }
#       
#       fullData <- hits
#     }
# 
#     fullData <- DT::datatable(fullData, extensions = 'Buttons',
#                                 escape = FALSE,
#                                 rownames = FALSE,
#                                 options = list(dom = 'Bfrtip',
#                                                buttons = 
#                                                  list('colvis', list(
#                                                    extend = 'collection',
#                                                    buttons = list(list(extend='csv',
#                                                                        filename = 'epHits'),
#                                                                   list(extend='excel',
#                                                                        filename = 'epHits'),
#                                                                   list(extend='pdf',
#                                                                        filename= 'epHits')),
#                                                    text = 'Download'
#                                                  )),
#                                                scrollX = TRUE,
#                                                pageLength = nrow(fullData),
#                                                order=list(list(2,'desc'))))
#     
#   })

})