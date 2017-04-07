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
  
###################################################################
   
###############################################################
# Graphs:   
  source("stackPlot.R",local=TRUE)$value
  source("boxPlot.R",local=TRUE)$value
  source("heatMap.R",local=TRUE)$value
################################################################


# # #############################################################    
#    
  output$mymap <- leaflet::renderLeaflet({

    map <- leaflet(height = "50px") %>%
      addProviderTiles("CartoDB.Positron") %>%
      setView(lng = -83.5, lat = 44.5, zoom=6)

    map

  })
  
  output$mapFooter <- renderUI({
    
    validate(
      need(!is.null(input$data), "Please select a data set")
    )
    
    chemicalSummary <- chemicalSummary()
    
    statsOfGroupOrdered <- chemicalSummary %>%
      group_by(site, date,category) %>%
      summarise(sumEAR = sum(EAR),
                hits = sum(EAR > hit_threshold)) %>%
      group_by(site,category) %>%
      summarise(mean = sum(mean(sumEAR) > hit_threshold),
                max = sum(max(sumEAR) > hit_threshold)) %>%
      data.frame()
    
    meanEARlogic <- input$meanEAR
    
    if(as.logical(meanEARlogic)){
      counts <- statsOfGroupOrdered$mean
    } else {
      counts <- statsOfGroupOrdered$max
    }
    
    if(input$radioMaxGroup == "1"){
      word <- "groups"
    } else if (input$radioMaxGroup == "2"){
      word <- "chemicals"
    } else {
      word <- "classes"
    }
    
    HTML(paste0("<h5>Size range represents number of ",word,
                " with hits. Ranges from ", min(counts,na.rm = TRUE)," - ", max(counts,na.rm = TRUE),"</h5>"))
    
  })
  
#   
#   observe({
#     
#     chemGroup <- chemicalSummary()
#     
#     # chemGroup <- filter(chemGroup, class == "Herbicide")
#     
#     meanEARlogic <- input$meanEAR
#     
#     maxEARWords <- ifelse(meanEARlogic,"meanEAR","maxEAR")
# 
#     if(input$radioMaxGroup == "1"){
#       typeWords <- "groups"
#     } else if (input$radioMaxGroup == "2"){
#       typeWords <- "chemicals"
#     } else {
#       typeWords <- "classes"
#     }
#       
#     statsOfGroupOrdered <- statsOfGroupOrdered()
#     
#     sumStat <- chemGroup %>%
#       group_by(site, date) %>%
#       summarise(sumEAR = sum(EAR),
#                 hits=as.numeric(any(hits > 0))) %>%
#       data.frame() %>%
#       group_by(site) %>%
#       summarise(nSamples = n(),
#                 meanEAR=ifelse(as.logical(meanEARlogic),mean(sumEAR,na.rm=TRUE),max(sumEAR,na.rm=TRUE)),
#                 freq=sum(hits)/n()) %>%
#       data.frame()
# 
#     siteToFind <- unique(sumStat$site)
# 
#     mapData <- right_join(stationINFO[,c("shortName", "Station.Name", "dec.lat.va","dec.long.va")], sumStat, by=c("shortName"="site"))
#     mapData <- left_join(mapData, statsOfGroupOrdered, by=c("shortName"="site", "nSamples"="nSamples"))
#     
#     mapData <- mapData[!is.na(mapData$dec.lat.va),]
#     mapData <- unique(mapData)
#     
#     col_types <- c("darkblue","dodgerblue","green4","gold1","orange","brown","red")
# 
#     if(as.logical(meanEARlogic)){
#       counts <- mapData$mean
#     } else {
#       counts <- mapData$max
#     }
#     
#     mapData$meanEAR[is.na(mapData$meanEAR)] <- 0
#     counts[is.na(counts)] <- 0
#     
#     if(length(siteToFind) > 1){
#       leg_vals <- unique(as.numeric(quantile(mapData$meanEAR, probs=c(0,0.01,0.1,0.25,0.5,0.75,0.9,.99,1), na.rm=TRUE)))
#       pal = colorBin(col_types, mapData$meanEAR, bins = leg_vals)
#       rad <-3*seq(1,4,length.out = 16)
#       if(sum(counts) == 0){
#         mapData$sizes <- rad[1]
#       } else {
#         mapData$sizes <- rad[as.numeric(cut(counts, breaks=16))]
#       }
#       
#     } else {
#       leg_vals <- unique(as.numeric(quantile(c(0,mapData$meanEAR), probs=c(0,0.01,0.1,0.25,0.5,0.75,0.9,.99,1), na.rm=TRUE)))
#       pal = colorBin(col_types, c(0,mapData$meanEAR), bins = leg_vals)
#       mapData$sizes <- 3
#     }
# 
#     # leg_vals <- c(0,0.0001,0.001,0.01,0.1,1,10)
#     # pal = colorBin(col_types, c(0,mapData$maxEAR), bins = leg_vals)
# 
#     map <- leafletProxy("mymap", data=mapData) %>%
#       clearMarkers() %>%
#       clearControls() %>%
#       addCircleMarkers(lat=~dec.lat.va, lng=~dec.long.va, 
#                  popup=paste0('<b>',mapData$Station.Name,"</b><br/><table>",
#                               "<tr><td>",maxEARWords,": </td><td>",sprintf("%.1f",mapData$meanEAR),'</td></tr>',
#                               "<tr><td>Number of Samples: </td><td>",mapData$nSamples,'</td></tr>',
#                               "<tr><td>Frequency: </td><td>",sprintf("%.1f",mapData$freq),'</td></tr>',
#                               "<tr><td>Number of ",typeWords," with hits: </td><td>",counts,'</td></tr>',
#                               '</table>') ,
#                  fillColor = ~pal(meanEAR), 
#                  # weight = 1,
#                  # color = "black",
#                  fillOpacity = 0.8, 
#                  radius = ~sizes, 
#                  stroke=FALSE,
#                  opacity = 0.8) 
#     
# 
#     if(length(siteToFind) > 1){
#       map <- addLegend(map,
#         position = 'bottomleft',
#         pal=pal,
#         values=~meanEAR,
#         opacity = 0.8,
#         labFormat = labelFormat(digits = 2), #transform = function(x) as.integer(x)),
#         title = ifelse(meanEARlogic,'Mean EAR','Max EAR'))
#         # title = "Herbicide Max EAR")
#     }
#     
#     map
#     
#   })
#   
# ############################################################# 
#   
  output$TableHeader <- renderUI({
    HTML(paste("<h4>", input$group,"-",input$data,": ",input$sites, "</h4>"))
  })



  output$BoxHeader <- renderUI({
    HTML(paste("<h4>", input$group,"-",input$data,": ",input$sites, "</h4>"))
  })

  output$TableHeaderColumns <- renderUI({
    HTML(paste("<h4>", input$data,": ",input$sites, "</h4>"))
  })

  output$siteHitText <- renderUI({

    if(input$sites == "All"){
      HTML(paste("<h4>Number of sites with hits</h4>"))
    } else {
      HTML(paste("<h4>Number of samples with hits</h4>"))
    }

  })
 
  output$nGroup <- renderUI({

    radio <- input$radioMaxGroup
    chemicalSummary <- chemicalSummary()
    siteToFind <- length(unique(chemicalSummary$site))

    if(radio == "2"){
      word <- "chemicals"
    } else if (radio == "3") {
      word <- "classes"
    } else if(radio == "4"){
      word <- "endPoints"
    } else {
      word <- "groups"
    }

    if(length(siteToFind) > 1){
      place <- "per site"
    } else {
      place <- ""
    }

    if(length(siteToFind) > 1){
      textUI <- paste("<h5>max = Maximum number of",word,"with hits per site</h5>",
                      "<h5>mean = Mean number of",word,"with hits per site</h5>",
                      "<h5>nSamples = Number of samples per site</h5>")
    } else {
      textUI <- paste("<h5>hits = Number of samples with hits </h5>",
                      "<h5>nSamples = Number of samples per site</h5>")
    }
    HTML(textUI)

  })
#   
#   output$endpointGraph <- renderPlot({ 
# 
#     filterBy <- input$epGroup
#     meanEARlogic <- as.logical(input$meanEAR)
#     hitThres <- hitThresValue()
#     
#     filterCat <- switch(as.character(input$radioMaxGroup),
#                         "1" = "choices",
#                         "2" = "chnm",
#                         "3" = "class")
#     
#     boxData <- chemicalSummary()
#   
#     graphData <- boxData %>%
#       # filter(!is.na(category)) %>%
#       group_by(site,date,category,endPoint) %>%
#       summarise(sumEAR=sum(EAR)) %>%
#       data.frame() %>%
#       group_by(site, category,endPoint) %>%
#       summarise(meanEAR=ifelse(meanEARlogic,mean(sumEAR),max(sumEAR))) %>%
#       data.frame() %>%
#       mutate(category=as.character(category)) 
# 
#     if(filterBy != "All"){
#       graphData <- graphData %>%
#         filter_(paste0("category == '", filterBy,"'"))
#       
#       countNonZero <- graphData %>%
#         group_by(endPoint) %>%
#         summarise(nonZero = as.character(sum(meanEAR>0)),
#                   hits = as.character(sum(meanEAR>hitThres)))
#       
#       countNonZero$hits[countNonZero$hits == "0"] <- ""
#       
#       namesToPlotEP <<- as.character(countNonZero$endPoint)
#       nSamplesEP <<- countNonZero$nonZero
#       nHitsEP <<- countNonZero$hits
#     }
#     
#     orderColsBy <- graphData %>%
#       group_by(endPoint) %>%
#       summarise(median = quantile(meanEAR[meanEAR != 0],0.5)) %>%
#       arrange(median)
#     
#     orderedLevelsEP <<- orderColsBy$endPoint
#     
#     if(any(is.na(orderColsBy$median))){
#       orderedLevelsEP <<- c(orderColsBy$endPoint[is.na(orderColsBy$median)],
#                           orderColsBy$endPoint[!is.na(orderColsBy$median)])
#     }
#     
#     graphData$endPoint <- factor(graphData$endPoint, levels = orderedLevelsEP)
#     
#     stackedPlot <- ggplot(graphData)+
#       scale_y_log10(paste(ifelse(meanEARlogic,"Mean","Maximum"), "EAR Per Site"),labels=fancyNumbers) +
#       geom_boxplot(aes(x=endPoint, y=meanEAR)) +
#       theme_minimal() +
#       xlab("") +
#       theme(axis.text.y = element_text(vjust = .25,hjust=1)) +
#       geom_hline(yintercept = hitThres, linetype="dashed", color="black")
#     
#     if(filterBy != "All"){
# 
#       ymin <<- 10^(ggplot_build(stackedPlot)$layout$panel_ranges[[1]]$y.range)[1]
#       ymax <<- 10^(ggplot_build(stackedPlot)$layout$panel_ranges[[1]]$y.range)[2]
# 
#       xmax <<- ggplot_build(stackedPlot)$layout$panel_ranges[[1]]$x.range[2]
#       xmin <<- ggplot_build(stackedPlot)$layout$panel_ranges[[1]]$x.range[1]
# 
#       stackedPlot <- stackedPlot +
#         geom_text(data=data.frame(), aes(x=namesToPlotEP, y=ymin,label=nSamplesEP),size=5) +
#         geom_text(data=data.frame(), aes(x=namesToPlotEP, y=ymax,label=nHitsEP),size=5)
# 
#       df1 <- data.frame(y = c(ymin,hitThres,ymax), text = c("# Non Zero","Hit Threshold","# Hits"), stringsAsFactors = FALSE)
# 
#       for(i in 1:3){
#         stackedPlot <- stackedPlot +
#           annotation_custom(
#             grob = textGrob(label = df1$text[i], gp = gpar(cex = 0.75)),
#             ymin = log10(df1$y[i]),      # Vertical position of the textGrob
#             ymax = log10(df1$y[i]),
#             xmin = xmax+0.05,  # Note: The grobs are positioned outside the plot area
#             xmax = xmax+0.05)
#       }
#     }
#     stackedPlot <- stackedPlot +
#       coord_flip()
#     
#     # Code to override clipping
#     stackedPlot <- ggplot_gtable(ggplot_build(stackedPlot))
#     stackedPlot$layout$clip[stackedPlot$layout$name == "panel"] <- "off"
#     
#     print(grid.draw(stackedPlot))
#     
#   })
#   
#   output$hitsTable <- DT::renderDataTable({    
#     
#     boxData <- chemicalSummary()
#     meanEARlogic <- as.logical(input$meanEAR)
#     hitThres <- hitThresValue()
#     
#     if(length(unique(boxData$site)) > 1){
#       tableData <- boxData %>%
#         group_by(site, choices, category, date) %>%
#         summarize(sumEAR = sum(EAR)) %>%
#         group_by(site, choices, category) %>%
#         summarize(meanEAR = ifelse(meanEARlogic, mean(sumEAR),max(sumEAR))) %>%
#           # hits = any(hits > 0)) %>% #is a hit when any EAR is greater than 0.1?
#         group_by(choices, category) %>%
#         summarize(nSites = sum(meanEAR>hitThres)) %>%
#         data.frame() 
#     } else {
#       tableData <- boxData %>%
#         group_by(choices, category, date)
#       
#       tableData <- tableData %>%
#         summarise(sumEAR=sum(EAR)) %>%
#         data.frame() %>%
#         group_by(choices, category) %>%
#         summarise(nSites = sum(sumEAR>hitThres))%>%
#         data.frame() 
#       
#     }
#     
#     if(input$radioMaxGroup != "1"){
#       tableData <- tableData %>%
#         reshape(idvar="choices",timevar="category", direction="wide") 
#       
#       names(tableData) <- gsub("nSites.","",names(tableData))
#       names(tableData)[1] <- "Groups"
#       
#       sumOfColumns <- colSums(tableData[-1],na.rm = TRUE)
#       orderData <- order(sumOfColumns,decreasing = TRUE) 
#       orderData <- orderData[sumOfColumns[orderData] != 0] + 1
#       
#       tableData <- tableData[,c(1,orderData)]
#       colors <- brewer.pal(9,"Blues") #"RdYlBu"
#       
#       groups <- tableData$Groups
#       
#       tableData <- tableData[!is.na(groups),-1,drop=FALSE]
#       rownames(tableData) <- groups[!is.na(groups)]
#       
#       cuts <- seq(0,max(as.matrix(tableData),na.rm=TRUE),length.out = 8)
#       
#       names(tableData)[names(tableData) == "Human Drug, Non Prescription"] <- "Human Drug"
#       names(tableData)[names(tableData) == "Flavor/Fragrance"] <- "Flavor / Fragrance"
#     } else {
#       tableData <- select(tableData, choices, nSites)
#       rownames(tableData) <- tableData$choices
#       tableData <- tableData[,-1,drop=FALSE]
#     }
#     
#     tableData1 <- DT::datatable(tableData, extensions = 'Buttons',
#                                 rownames = TRUE,
#                                 options = list(scrollX = TRUE,
#                                                dom = 'Bfrtip',
#                                                buttons = 
#                                                  list('colvis', list(
#                                                    extend = 'collection',
#                                                    buttons = list(list(extend='csv',
#                                                                        filename = 'siteHits'),
#                                                                   list(extend='excel',
#                                                                        filename = 'siteHits'),
#                                                                   list(extend='pdf',
#                                                                        filename= 'siteHits')),
#                                                    text = 'Download',
#                                                    filename= 'test'
#                                                  )),
#                                  pageLength = nrow(tableData),
#                                  order=list(list(1,'desc'))))
#     if(input$radioMaxGroup != "1"){
#       for(i in 1:ncol(tableData)){
#         tableData1 <- formatStyle(tableData1, columns = names(tableData)[i], 
#                     backgroundColor = styleInterval(cuts = cuts,values = colors),
#                     color = styleInterval(0.75*max(tableData,na.rm=TRUE),values = c("black","white")),
#                     `font-size` = '17px')
#       }
#     }
#     tableData1
#     
#   })
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