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

endPointInfo <- endPointInfo

endPointInfo <- endPointInfo[!(endPointInfo$assay_source_name == "ATG" & endPointInfo$signal_direction == "loss"),]
endPointInfo <- endPointInfo[!(endPointInfo$assay_source_name == "NVS" & endPointInfo$signal_direction == "gain"),]
endPointInfo <- endPointInfo[endPointInfo$assay_component_endpoint_name != "TOX21_p53_BLA_p3_ratio",]
endPointInfo <- endPointInfo[endPointInfo$assay_component_endpoint_name != "TOX21_p53_BLA_p2_viability",]

choicesPerGroup <- apply(endPointInfo, 2, function(x) length(unique(x[!is.na(x)])))

choicesPerGroup <- which(choicesPerGroup > 6 & choicesPerGroup < 100)

endPointInfo <- endPointInfo[,c(39,as.integer(choicesPerGroup))]

cleanUpNames <- endPointInfo$intended_target_family
cleanUpNames <- stri_trans_totitle(cleanUpNames)
cleanUpNames[grep("Dna",cleanUpNames)] <- "DNA Binding"
cleanUpNames[grep("Cyp",cleanUpNames)] <- "CYP"
cleanUpNames[grep("Gpcr",cleanUpNames)] <- "GPCR"
endPointInfo$intended_target_family <- cleanUpNames

sitesOrdered <- c("StLouis","Pigeon","Nemadji","WhiteWI","Bad","Montreal","PresqueIsle",
                  "Ontonagon","Sturgeon","Tahquamenon","Manistique","Escanaba","Ford","Cheboygan2","Indian",
                  "Menominee","Peshtigo","Oconto","Fox","Manistee","Manitowoc","PereMarquette","Sheboygan",
                  "WhiteMI","Muskegon","MilwaukeeMouth","GrandMI","Kalamazoo2","PawPaw",
                  "StJoseph","IndianaHC","Burns","ThunderBay","AuSable","Rifle",
                  "Saginaw","BlackMI","Clinton","Rouge","HuronMI","Raisin","Maumee",
                  "Portage","Sandusky","HuronOH","Vermilion","BlackOH","Rocky","Cuyahoga","GrandOH",
                  "Cattaraugus","Tonawanda","Genesee","Oswego","BlackNY","Oswegatchie","Grass","Raquette","StRegis")

pathToApp <- system.file("extdata", package="toxEval")

stationINFO <- readRDS(file.path(pathToApp,"sitesOWC.rds"))
summaryFile <- readRDS(file.path(pathToApp,"summary.rds"))
endPoint <- readRDS(file.path(pathToApp,"endPoint.rds"))

df2016 <- readRDS(file.path(pathToApp,"df2016.rds"))

choicesPerGroup <- apply(endPointInfo[,-1], 2, function(x) length(unique(x)))
groupChoices <- paste0(names(choicesPerGroup)," (",choicesPerGroup,")")

interl <- function (a,b) {
  n <- min(length(a),length(b))
  p1 <- as.vector(rbind(a[1:n],b[1:n]))
  p2 <- c(a[-(1:n)],b[-(1:n)])
  c(p1,p2)
}

fancyNumbers <- function(n){
  nNoNA <- n[!is.na(n)]
  x <-gsub(pattern = "1e",replacement = "10^",x = format(nNoNA, scientific = TRUE))
  exponents <- as.numeric(sapply(strsplit(x, "\\^"), function(j) j[2]))
  base <- ifelse(exponents == 0, "1", ifelse(exponents == 1, "10","10^"))
  exponents[exponents == 0 | exponents == 1] <- ""
  textNums <- rep(NA, length(n))  
  textNums[!is.na(n)] <- paste0(base,exponents)

  textReturn <- parse(text=textNums)
  return(textReturn)
}

fancyNumbers2 <- function(n){
  textReturn <-  signif(n,digits = 2)
  textReturn <- as.character(textReturn)
  textReturn[length(textReturn)] <- paste(">",textReturn[length(textReturn)])
  textReturn[1] <- paste("<",textReturn[1])
  return(textReturn)
}

createLink <- function(cas,ep, hits) {
  paste0('<a href="http://actor.epa.gov/dashboard/#selected/',cas,"+",ep,'" target="_blank" >',hits,'</a>')
}

shinyServer(function(input, output,session) {
  
  choices <- reactive({
    if (input$data == "Duluth" ){
      duluthSites <- readRDS(file.path(pathToApp,"sitesOWC.rds"))
      choices =  c("All",duluthSites$shortName)
    } else if(input$data == "NPS"){
      npsSites <- readRDS(file.path(pathToApp,"npsSite.rds"))
      choices =  c("All",npsSites$shortName)      
    } else {
      choices =  c("All","Potential 2016",summaryFile$site)
    }
    
    choices
  })

  observe({
    updateSelectInput(session, "sites", choices = choices())
  })

  assayDF <- eventReactive(input$pickAssay, ignoreNULL = FALSE, {
    
    ep <- groupDF()
    assays <- input$assay
    
    endPointInfoSub <- distinct(ep) %>%
      rename(endPoint=assay_component_endpoint_name) %>%
      data.frame() 
    
    index <- unique(grep(paste(assays,collapse="|"), 
                         endPointInfoSub$endPoint))
    
    endPointInfoSub <- endPointInfoSub[index,]
    
    endPointInfoSub
    
  })
  
  groupDF <- eventReactive(input$changeAnn, ignoreNULL = FALSE, {
    groupCol <- input$groupCol

    ep <- data.frame(endPointInfo[,c("assay_component_endpoint_name", groupCol)])
    ep
  })
  
  hitThresValue <- eventReactive(input$changeHit, ignoreNULL = FALSE, {
    hitThresValue <- input$hitThres
    hitThresValue
  })

  observe({
    
    ep <- groupDF()
    ep <- ep[!is.na(ep[,2]),]
    ep <- ep[ep[,2] != "NA",]
    
    orderBy <- ep[,2]
    orderNames <- names(table(orderBy))
    nEndPoints <- as.integer(table(orderBy))
    
    df <- data.frame(orderNames,nEndPoints,stringsAsFactors = FALSE) %>%
      arrange(desc(nEndPoints))

    dropDownHeader <- c(paste0(df$orderNames," (",df$nEndPoints,")"))

    selChoices <- df$orderNames

    if(names(ep)[2] == "intended_target_family"){
      selChoices <- selChoices[!(selChoices %in% c("Cell Cycle","Background Measurement","Cell Morphology"))]
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
      valueText <- c("All",unique(chemicalSummary$class))
    } else if(input$radioMaxGroup == 1){
      valueText <- c("All",unique(chemicalSummary$choices))
    }
    
    updateSelectInput(session, "epGroup", choices = valueText)
  })
  
#############################################################   
  chemicalSummary <- reactive({
    
    ep <- assayDF()
    path <- pathToApp
    
    groupCol <- names(ep)[2]
    
    radioMaxGroup <- input$radioMaxGroup
    group <- input$group
    
    if(!any(group %in% names(table(ep[,2])))){
      group <- names(table(ep[,2]))
    }
    
    # validate(
    #   need(length(input$assay) > 0, 'Check at least one assay')
    # )
    
    if (input$data == "Water Sample"){
      chemicalSummary <- readRDS(file.path(path,"chemicalSummaryV2.rds"))
      # chemicalSummary <- readRDS(file.path(path,"chemicalSummary_noFlags.rds"))
      stationINFO <<- readRDS(file.path(path,"sitesOWC.rds"))
    } else if (input$data == "Passive Samples"){
      chemicalSummary <- readRDS(file.path(path,"chemicalSummaryPassive.rds"))
      chemicalSummary$date <- as.POSIXct(as.Date(paste0(chemicalSummary$year,"-01-01")))
      if(input$year == "2014"){
        chemicalSummary <- filter(chemicalSummary, date > as.POSIXct(as.Date("2011-01-01")))
      } else if (input$year == "2010"){
        chemicalSummary <- filter(chemicalSummary, date < as.POSIXct(as.Date("2011-01-01")))
      }
      
      stationINFO <<- readRDS(file.path(path,"sitesOWC.rds"))
    } else if (input$data == "Duluth"){
      chemicalSummary <- readRDS(file.path(path,"chemSummeryDL.rds"))
      stationINFO <<- readRDS(file.path(path,"sitesDuluth.rds"))
    } else if (input$data == "NPS"){
      chemicalSummary <- readRDS(file.path(path,"chemNPS.rds"))
      stationINFO <<- readRDS(file.path(path,"npsSite.rds"))
    } else if (input$data == "Detection Limits"){
      chemicalSummary <- readRDS(file.path(path,"chemicalSummaryDetectionLevels.rds"))
      stationINFO <<- readRDS(file.path(path,"sitesOWC.rds"))
    }
    
    chemicalSummary <- chemicalSummary %>%
      filter(endPoint %in% ep$endPoint) %>%
      data.table() %>%
      left_join(data.table(ep), by="endPoint") %>%
      data.frame() %>%
      select_("hits","EAR","chnm","class","date","choices"=groupCol,"site","endPoint","endPointValue","casrn") %>%
      left_join(stationINFO[,c("fullSiteID","shortName")], by=c("site"="fullSiteID")) %>%
      select(-site) %>%
      rename(site=shortName)
    
    names(ep)[2] <- "groupC"
    
    endPointInfoSub <- filter(ep, groupC %in%  group)
    
    chemicalSummary <- filter(chemicalSummary, endPoint %in% endPointInfoSub$endPoint)

  
    if(radioMaxGroup == "2"){
      chemicalSummary <- chemicalSummary %>%
        mutate(category=chnm)
    } else if(radioMaxGroup == "3"){
      chemicalSummary <- chemicalSummary %>%
        mutate(category=class)
    } else if(radioMaxGroup == "4"){
      chemicalSummary <- chemicalSummary %>%
        mutate(category=endPoint) 
    } else {
      chemicalSummary <- chemicalSummary %>%
        mutate(category = choices)         
    }
    
    if (input$sites == "Potential 2016"){
      chemicalSummary <-  chemicalSummary %>%
        filter(site %in% df2016$Site)
    } else if (input$sites != "All"){
      if(any(input$sites %in% chemicalSummary$site)){
        chemicalSummary <-  chemicalSummary %>%
          filter(site %in% input$sites)
      }
    }
    
    if(length(unique(chemicalSummary$site)) == 1){
      chemicalSummary <- chemicalSummary %>%
        mutate(date = factor(as.numeric(date) - min(as.numeric(date))))
    }
    
    chemicalSummary 
    
  })
  
  statsOfColumn <- reactive({
    
    chemicalSummary <- chemicalSummary()
    meanEARlogic <- input$meanEAR

    radio <- input$radioMaxGroup
    statsOfColumn <- chemicalSummary 

    siteToFind <- unique(statsOfColumn$site)
    
    if(length(siteToFind) == 1){

      if(radio == "2"){
        statsOfColumn <- mutate(statsOfColumn, site = chnm) 
      } else if (radio == "3"){
        statsOfColumn <- mutate(statsOfColumn, site = class) 
      } else {
        statsOfColumn <- mutate(statsOfColumn, site = choices) 
      }
    } 

    statsOfColumn <- statsOfColumn %>%
      group_by(site, date, category) %>%
      summarise(sumEAR = sum(EAR),
                nHits = sum(hits)) %>%
      group_by(site, category) %>%
      summarise(maxEAR = ifelse(meanEARlogic, mean(sumEAR), max(sumEAR)),
                freq = sum(nHits > 0)/n()) %>%
      data.frame()
    
    if(!(length(siteToFind) == 1)){ #  & radio == "1"
      statsOfColumn <- statsOfColumn %>%
        gather(calc, value, -site, -category) %>%
        unite(choice_calc, category, calc, sep=" ") %>%
        spread(choice_calc, value)        
    }

    statsOfColumn
    
  })
  
  statsOfGroupOrdered <- reactive({

    statsOfGroup <- chemicalSummary()   
    hitThres <- hitThresValue()
    
    siteToFind <- unique(statsOfGroup$site)
    
    statsOfGroupOrdered <- statsOfGroup %>%
      group_by(site, date,category) %>%
      summarise(sumEAR = sum(EAR)) %>%
      group_by(site,category) %>%
      summarise(max = sum(sumEAR > hitThres),
                mean = sum(mean(sumEAR) > hitThres),
                hit = as.numeric(any(sumEAR > hitThres)),
                nSamples = n())%>%
      data.frame()

    if(length(siteToFind) > 1){
      statsOfGroupOrdered <- statsOfGroupOrdered %>%
        group_by(site) %>%
        summarise(max=sum(max > 0),
                  mean=sum(mean > 0),
                  # total=sum(hit>=1),
                  nSamples = median(nSamples))

    }

    statsOfGroupOrdered
      
  }) 
  
# ############################################################# 

  output$tableGroupSumm <- DT::renderDataTable({        

    statsOfGroupOrdered <- statsOfGroupOrdered()

    siteToFind <- unique(statsOfGroupOrdered$site)

    if(length(siteToFind) > 1){
      
      colToSort <- 1
      
      meanChem <- grep("mean",names(statsOfGroupOrdered))
      maxChem <- grep("max",names(statsOfGroupOrdered))

      tableGroup <- DT::datatable(statsOfGroupOrdered,  extensions = 'Buttons',
                                    rownames = FALSE,
                                    options = list(
                                                  pageLength = nrow(statsOfGroupOrdered), 
                                                  order=list(list(colToSort,'desc')),
                                                 dom = 'Bfrtip',
                                                 buttons =
                                                   list('colvis', list(
                                                     extend = 'collection',
                                                     buttons = list(list(extend='csv',
                                                                         filename = 'hitCount'),
                                                                    list(extend='excel',
                                                                         filename = 'hitCount'),
                                                                    list(extend='pdf',
                                                                         filename= 'hitCount')),
                                                     text = 'Download')
                                                   )
                                                 ))
                                    # options = list(dom = 'ft',


      tableGroup <- formatStyle(tableGroup, names(statsOfGroupOrdered)[maxChem], 
                                background = styleColorBar(range(statsOfGroupOrdered[,maxChem],na.rm=TRUE), 'goldenrod'),
                                backgroundSize = '100% 90%',
                                backgroundRepeat = 'no-repeat',
                                backgroundPosition = 'center' ) 
      
      tableGroup <- formatStyle(tableGroup, names(statsOfGroupOrdered)[meanChem], 
                                background = styleColorBar(range(statsOfGroupOrdered[,meanChem],na.rm=TRUE), 'wheat'),
                                backgroundSize = '100% 90%',
                                backgroundRepeat = 'no-repeat',
                                backgroundPosition = 'center') 

      } else {
        
        tableGroup <- DT::datatable(statsOfGroupOrdered[,c("category","max","nSamples")], extensions = 'Buttons', 
                                    colnames = c('hits' = 2),
                                    rownames = FALSE,
                                    options = list(
                                      pageLength = nrow(statsOfGroupOrdered), 
                                      order=list(list(1,'desc')),
                                      dom = 'Bfrtip',
                                      buttons =
                                        list('colvis', list(
                                          extend = 'collection',
                                          buttons = list(list(extend='csv',
                                                              filename = 'hitCount'),
                                                         list(extend='excel',
                                                              filename = 'hitCount'),
                                                         list(extend='pdf',
                                                              filename= 'hitCount')),
                                          text = 'Download')
                                        )
                                    ))

        tableGroup <- formatStyle(tableGroup, "hits", 
                                  background = styleColorBar(range(statsOfGroupOrdered[,"max"],na.rm=TRUE), 'goldenrod'),
                                  backgroundSize = '100% 90%',
                                  backgroundRepeat = 'no-repeat',
                                  backgroundPosition = 'center' ) 

      }
      
      tableGroup
    
  })
  
  output$tableSumm <- DT::renderDataTable({
    
    statCol <- statsOfColumn()
    meanEARlogic <- as.logical(input$meanEAR)
    
    colToSort <- 1
    if("nSamples" %in% names(statCol)){
      colToSort <- 2
    }
    freqCol <- grep("freq",names(statCol))
    maxEARS <- grep("maxEAR",names(statCol))
    
    ignoreIndex <- which(names(statCol) %in% c("site","nSamples"))
    
    statCol <- statCol[,c(ignoreIndex,c(maxEARS,freqCol)[order(c(maxEARS,freqCol))])]
    
    maxEARS <- grep("maxEAR",names(statCol))
    
    MaxEARSordered <- order(apply(statCol[,maxEARS, drop = FALSE], 2, max),decreasing = TRUE)
    
    if(length(maxEARS) > 9){
      statCol <- statCol[,c(ignoreIndex,interl(maxEARS[MaxEARSordered[1:9]],(maxEARS[MaxEARSordered[1:9]]-1)))]
      maxEARS <- maxEARS[1:9]
    } else {
      if(length(maxEARS) != 1){
        statCol <- statCol[,c(ignoreIndex,interl(maxEARS[MaxEARSordered],(maxEARS[MaxEARSordered]-1)))]
      }
    }
    
    freqCol <- grep("freq",names(statCol))
    maxEARS <- grep("maxEAR",names(statCol))
    
    if(meanEARlogic){
      names(statCol)[maxEARS] <- gsub("max","mean",names(statCol)[maxEARS])
    }
    
    
    colors <- brewer.pal(length(maxEARS),"Blues") #"RdYlBu"
    tableSumm <- DT::datatable(statCol, extensions = 'Buttons', 
                               rownames = FALSE,
                               options = list(#dom = 'ft',
                                              dom = 'Bfrtip',
                                              buttons =
                                                list('colvis', list(
                                                  extend = 'collection',
                                                  buttons = list(list(extend='csv',
                                                                      filename = 'hitStats'),
                                                                 list(extend='excel',
                                                                      filename = 'hitStats'),
                                                                 list(extend='pdf',
                                                                      filename= 'hitStats')),
                                                  text = 'Download'
                                                  )
                                                ),
                                              scrollX = TRUE,
                                              pageLength = nrow(statCol),
                                              order=list(list(colToSort,'desc'))))
    
    tableSumm <- formatRound(tableSumm, names(statCol)[-ignoreIndex], 2) 
    
    for(i in 1:length(maxEARS)){
      tableSumm <- formatStyle(tableSumm, 
                               names(statCol)[maxEARS[i]], 
                               backgroundColor = colors[i])
      tableSumm <- formatStyle(tableSumm, 
                               names(statCol)[freqCol[i]], 
                               backgroundColor = colors[i])
      
      tableSumm <- formatStyle(tableSumm, names(statCol)[maxEARS[i]], 
                               background = styleColorBar(range(statCol[,names(statCol)[maxEARS[i]]],na.rm = TRUE), 'goldenrod'),
                               backgroundSize = '100% 90%',
                               backgroundRepeat = 'no-repeat',
                               backgroundPosition = 'center' ) 
      tableSumm <- formatStyle(tableSumm, names(statCol)[freqCol[i]], 
                               background = styleColorBar(range(statCol[,names(statCol)[freqCol[i]]],na.rm = TRUE), 'wheat'),
                               backgroundSize = '100% 90%',
                               backgroundRepeat = 'no-repeat',
                               backgroundPosition = 'center') 
      
    }
    
    tableSumm
    
  })
  
# #############################################################
   
  output$stackBarGroup <- renderPlot({ 
    
    graphData <- graphData()
    meanEARlogic <- as.logical(input$meanEAR)
    catType = as.numeric(input$radioMaxGroup)
    siteToFind <- unique(graphData$site)
    
    cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    cbValues <- colorRampPalette(cbPalette)(length(levels(graphData$reorderedCat)))
    set.seed(4)
    cbValues <- sample(cbValues)
    
    siteLimits <- stationINFO %>%
      filter(shortName %in% unique(graphData$site))
    
    if(length(siteToFind) > 1){
      
      if(all(siteLimits$shortName %in% sitesOrdered)){
        siteLimits <- mutate(siteLimits, shortName = factor(shortName, levels=sitesOrdered[sitesOrdered %in% siteLimits$shortName]))
      } else {
        siteLimits <- mutate(siteLimits, shortName = factor(shortName))
      }
      
      if(catType != 2){
        upperPlot <- ggplot(graphData, aes(x=site, y=meanEAR, fill = reorderedCat))
      } else {
        upperPlot <- ggplot(graphData, aes(x=site, y=meanEAR, fill = class))
      }
      
      upperPlot <- upperPlot +
        geom_bar(stat="identity") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.25,
                                         colour=siteLimits$lakeColor)) +
        scale_x_discrete(limits=levels(siteLimits$shortName),drop=FALSE) +
        xlab("") +
        ylab(paste(ifelse(meanEARlogic,"Mean","Maximum"), "EAR Per Site")) +
        scale_fill_manual(values = cbValues, drop=FALSE) + 
        guides(fill=FALSE) 

      
    } else {
      
      chemGroupBPOneSite <- chemicalSummary() %>%
        select(-site) 
      
      chemGroupBPOneSite <- mutate(chemGroupBPOneSite, category = factor(category,levels=orderedLevels))
      
      upperPlot <- ggplot(chemGroupBPOneSite, aes(x=date, y=EAR, fill = category)) +
        geom_bar(stat="identity")+
        theme_minimal() +
        theme(axis.text.x=element_blank(),
              axis.ticks=element_blank())+
        xlab("Individual Samples") + 
        ylab("EAR") +
        scale_fill_discrete("", drop=FALSE) +
        scale_fill_manual(values = cbValues, drop=FALSE) +
        guides(fill=FALSE)

    }
    ggsave("stackPlot.png",upperPlot,bg = "transparent")
    
    print(upperPlot)
  })

  output$graphGroup <- renderPlot({ 
    
    hitThres <- hitThresValue()
    graphData <- graphData()
    meanEARlogic <- as.logical(input$meanEAR)
    catType = as.numeric(input$radioMaxGroup)
    siteToFind <- unique(graphData$site)
    boxPlotLogic <-  !(input$data == "Passive Samples" & length(siteToFind) == 1) 
    
    countNonZero <- graphData %>%
      group_by(category) %>%
      summarise(nonZero = as.character(sum(meanEAR>0)),
                hits = as.character(sum(meanEAR>hitThres)))
    
    countNonZero$hits[countNonZero$hits == "0"] <- ""
    
    namesToPlot <<- as.character(countNonZero$category)
    nSamples <<- countNonZero$nonZero
    nHits <<- countNonZero$hits
    
    if(length(siteToFind) > 1){
      
      lowerPlot <- ggplot(graphData)+
        scale_y_log10(paste(ifelse(meanEARlogic,"Mean","Maximum"), "EAR Per Site"),labels=fancyNumbers) # fance labels
      labelsText <<- "nSites"
      
      if(!boxPlotLogic){
        lowerPlot <- lowerPlot + geom_point(aes(x=reorderedCat, y=meanEAR, color=reorderedCat, size=3))
      } else {
        if(catType  != 2){
          lowerPlot <- lowerPlot + 
            geom_boxplot(aes(x=reorderedCat, y=meanEAR, fill=reorderedCat),lwd=0.1,outlier.size=1) 
          
        } else {
          lowerPlot <- lowerPlot + 
            geom_boxplot(aes(x=reorderedCat, y=meanEAR, fill=class),lwd=0.1,outlier.size=1) 
        }
      }
      
    } else {
      
      lowerPlot <- ggplot(graphData)+
        scale_y_log10("Sum of EAR",labels=fancyNumbers) # labels=fancyNumbers
      
      labelsText <<- "nSamples"
      
      if(!boxPlotLogic){
        lowerPlot <- lowerPlot + geom_point(aes(x=reorderedCat, y=meanEAR, color=reorderedCat, size=3))
      } else {
        lowerPlot <- lowerPlot + 
          geom_boxplot(aes(x=reorderedCat, y=meanEAR, fill=reorderedCat),lwd=0.1,outlier.size=1) 
      }
    }
    
    text.size <<- ifelse(length(orderedLevels) > 15, 13, 16)
    
    cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    cbValues <- colorRampPalette(cbPalette)(length(levels(graphData$reorderedCat)))
    set.seed(4)
    cbValues <- sample(cbValues)
    
    lowerPlot <- lowerPlot + 
      theme_minimal() +
      theme(axis.text.y = element_text(size=text.size, color = "black", margin = margin(0,-10,0,0),vjust = 0.2), 
            axis.text.x = element_text(size=16, color = "black", vjust = 0, margin = margin(-15,0,0,0)),
            axis.title = element_text(size=16)) +
      xlab("") +
      geom_hline(yintercept = hitThres, linetype="dashed", color="black",lwd=0.25)  +
      scale_x_discrete(drop=FALSE) +
      scale_fill_discrete(drop=FALSE)
    
    if(catType  != 2 | length(siteToFind) == 1){
      lowerPlot <- lowerPlot + 
        theme(legend.position = "none")
      
    } else {
      lowerPlot <- lowerPlot + 
        theme(legend.key = element_blank(),
              legend.title = element_blank(),
              # legend.text = element_text(size=9),
              legend.justification = c(1, 0), 
              legend.position = c(1, 0),
              legend.background = element_rect(colour = 'black', fill = 'white')) 
    }
    
    ymin <<- 10^(ggplot_build(lowerPlot)$panel$ranges[[1]]$y.range)[1]
    ymax <<- 10^(ggplot_build(lowerPlot)$panel$ranges[[1]]$y.range)[2]
    
    xmax <<- ggplot_build(lowerPlot)$panel$ranges[[1]]$x.range[2]
    xmin <<- ggplot_build(lowerPlot)$panel$ranges[[1]]$x.range[1]
    
    lowerPlot <- lowerPlot + 
      geom_text(data=data.frame(), aes(x=namesToPlot, y=ymin,label=nSamples),size=5)  +
      geom_text(data=data.frame(), aes(x=namesToPlot, y=ymax,label=nHits),size=5) 
    
    df1 <- data.frame(y = c(ymin,hitThres,ymax), text = c("# Non Zero","Hit Threshold","# Hits"), stringsAsFactors = FALSE)
    
    for(i in 1:3){
      lowerPlot <- lowerPlot + 
        annotation_custom(
          grob = textGrob(label = df1$text[i], gp = gpar(cex = 0.75)),
          ymin = log10(df1$y[i]),      # Vertical position of the textGrob
          ymax = log10(df1$y[i]),
          xmin = xmax+0.05,  # Note: The grobs are positioned outside the plot area
          xmax = xmax+0.05)
    }
    
    lowerPlot <- lowerPlot + 
      scale_fill_manual(values = cbValues, drop=FALSE) +
      coord_flip() 
    
    if(catType == 2 & length(siteToFind) > 1){
      lowerPlot <- lowerPlot +
        scale_color_manual(values = cbValues) 
    }
    
    # Code to override clipping
    lowerPlot <- ggplot_gtable(ggplot_build(lowerPlot))
    lowerPlot$layout$clip[lowerPlot$layout$name == "panel"] <- "off"

    ggsave("boxPlot.png",lowerPlot,bg = "transparent")
    
    print(grid.draw(lowerPlot))
  })
  
  output$downloadBoxPlot <- downloadHandler(
    
    filename = function() {
      "boxPlot.png"
    },
    content = function(file) {
      file.copy("boxPlot.png", file)
    }
  )
  
  output$downloadStackPlot <- downloadHandler(
    
    filename = function() {
      "stackPlot.png"
    },
    content = function(file) {
      file.copy("stackPlot.png", file)
    }
  )
  
  output$downloadHeatPlot <- downloadHandler(
    
    filename = function() {
      "heatPlot.png"
    },
    content = function(file) {
      file.copy("heatPlot.png", file)
    }
  )
  
  output$graphHeat <- renderPlot({ 
    
    graphData <- graphData()
    meanEARlogic <- as.logical(input$meanEAR)
    catType = as.numeric(input$radioMaxGroup)
    
    siteToFind <- unique(graphData$site)
    
    siteLimits <- stationINFO %>%
      filter(shortName %in% unique(graphData$site))
    
    if(all(siteLimits$shortName %in% sitesOrdered)){
      siteLimits <- mutate(siteLimits, shortName = factor(shortName, levels=sitesOrdered[sitesOrdered %in% siteLimits$shortName]))
    } else {
      siteLimits <- mutate(siteLimits, shortName = factor(shortName))
    }
    
    if(catType == 2 & length(siteToFind) > 1){

      levels(graphData$class)[levels(graphData$class) == "Human Drug, Non Prescription"] <- "Human Drug"
      levels(graphData$class)[levels(graphData$class) == "Antimicrobial Disinfectant"] <- "Antimicrobial"
      levels(graphData$class)[levels(graphData$class) == "Detergent Metabolites"] <- "Detergent"
      
      orderClass <- graphData %>%
        group_by(class,category) %>%
        summarise(median = median(meanEAR[meanEAR != 0])) %>%
        data.frame() %>%
        arrange(desc(median)) %>%
        filter(!duplicated(class)) %>%
        arrange(median) 
      
      orderChem <- graphData %>%
        group_by(category,class) %>%
        summarise(median = quantile(meanEAR[meanEAR != 0],0.5)) %>%
        data.frame() %>%
        mutate(class = factor(class, levels=orderClass$class)) %>%
        arrange(class, median)
      
      orderClass <- mutate(orderClass, class = factor(class, levels=levels(orderChem$class)))
      
      orderedLevels <- orderChem$category[!is.na(orderChem$median)] 
      
      graphData$class <- factor(graphData$class, levels=levels(orderClass$class))
      graphData$reorderedCat <- factor(graphData$category, levels=orderedLevels)
      
      heat <- ggplot(data = graphData) +
        geom_tile(aes(x = site, y=reorderedCat, fill=meanEAR)) +
        theme(axis.text.x = element_text(colour=siteLimits$lakeColor,
                                         angle = 90,vjust=0.5,hjust = 1)) +
        scale_x_discrete(limits=levels(siteLimits$shortName),drop=FALSE) +
        ylab("") +
        xlab("") +
        labs(fill=paste(ifelse(meanEARlogic, "Mean","Maximum")," EAR")) +
        scale_fill_gradient( guide = "legend",
                             trans = 'log',
                             low = "white", high = "steelblue",
                             breaks=c(0.00001,0.0001,0.001,0.01,0.1,1,5),
                             na.value = 'lightgrey',labels=fancyNumbers2) +
        facet_grid(class ~ .,scales="free_y", space="free_y") +
        theme(strip.text.y = element_text(angle=0, hjust=0), 
              strip.background = element_rect(fill="white"),
              panel.margin.y=unit(0.05, "lines"))
      
    } else {
      
      graphData$reorderedCat <- factor(graphData$category, levels=rev(levels(graphData$reorderedCat)))
      
      heat <- ggplot(data = graphData) +
        geom_tile(aes(x = site, y=reorderedCat, fill=meanEAR)) +
        theme(axis.text.x = element_text(colour=siteLimits$lakeColor,
                                         angle = 90,vjust=0.5,hjust = 1)) +
        scale_x_discrete(limits=levels(siteLimits$shortName),drop=FALSE) +
        ylab("") +
        xlab("") +
        labs(fill=paste(ifelse(meanEARlogic, "Mean","Maximum")," EAR")) +
        scale_fill_gradient( guide = "legend",
                             trans = 'log',
                             low = "white", high = "steelblue",
                             breaks=c(0.00001,0.0001,0.001,0.01,0.1,1,5),
                             na.value = 'lightgrey',labels=fancyNumbers2)
    }

    ggsave("heatPlot.png",heat,bg = "transparent")
    
    print(heat)
  })
  
  output$graphHeat.ui <- renderUI({
    heightOfGraph <- 500
    if(as.numeric(input$radioMaxGroup) == 2){
      heightOfGraph <- 800
    }
    plotOutput("graphHeat", height = heightOfGraph)
  })
  
  output$graphGroup.ui <- renderUI({
    heightOfGraph <- 500
    if(as.numeric(input$radioMaxGroup) == 2){
      heightOfGraph <- 800
    }
    plotOutput("graphGroup", height = heightOfGraph)
  })
  
  graphData <- reactive({
    
    ep <- assayDF() 
    columnName <- names(ep)[2]
    
    meanEARlogic <- as.logical(input$meanEAR)
    catType <- as.numeric(input$radioMaxGroup)

    boxData <- chemicalSummary()
    siteToFind <- unique(boxData$site)
    
    if(catType == 1){
      orderNames <- names(table(ep[,2]))
    }
      
    if(catType == 2 & length(siteToFind) > 1){
        
        graphData <- chemicalSummary() %>%
          group_by(site,date,category,class) %>%
          summarise(sumEAR=sum(EAR)) %>%
          data.frame() %>%
          group_by(site, category,class) %>%
          summarise(meanEAR=ifelse(meanEARlogic,mean(sumEAR),max(sumEAR))) %>%
          data.frame() %>%
          mutate(category=as.character(category))
        
        orderClass <- graphData %>%
          group_by(class,category) %>%
          summarise(median = median(meanEAR[meanEAR != 0])) %>%
          data.frame() %>%
          arrange(desc(median)) %>%
          filter(!duplicated(class)) %>%
          arrange(median)
        
        orderChem <- graphData %>%
          group_by(category,class) %>%
          summarise(median = quantile(meanEAR[meanEAR != 0],0.5)) %>%
          data.frame() %>%
          mutate(class = factor(class, levels=orderClass$class)) %>%
          arrange(class, median)
        
        orderedLevels <<- orderChem$category # [!is.na(orderChem$median)] #The !is.na takes out any category that was all censo
        
        if(any(is.na(orderChem$median))){
          orderedLevels <<- c(orderChem$category[is.na(orderChem$median)],
                              orderChem$category[!is.na(orderChem$median)])
        }
        
        
        graphData <- mutate(graphData, reorderedCat = factor(as.character(category), levels=orderedLevels)) %>%
          mutate(class = factor(class, levels=rev(orderClass$class)))
      
      } else {
        
        graphData <- boxData %>%
          group_by(site,date,category) %>%
          summarise(sumEAR=sum(EAR)) %>%
          data.frame() 
        
        if (length(siteToFind) > 1){
    
          graphData <- graphData %>%
            group_by(site, category) %>%
            summarise(meanEAR=ifelse(meanEARlogic,mean(sumEAR),max(sumEAR))) %>%
            data.frame() %>%
            mutate(category=as.character(category)) 
          
        } else {
        
          graphData <- graphData %>%
            rename(meanEAR=sumEAR) %>%
            mutate(category=as.character(category)) 
          
        }
        
        orderColsBy <- graphData %>%
          mutate(category = as.character(category)) %>%
          group_by(category) %>%
          summarise(median = median(meanEAR[meanEAR != 0])) %>%
          arrange(median)
        
        orderedLevels <<- orderColsBy$category
        
        if(any(is.na(orderColsBy$median))){
          orderedLevels <<- c(orderColsBy$category[is.na(orderColsBy$median)],
                              orderColsBy$category[!is.na(orderColsBy$median)])
        }
        
        if(catType == 1){
          orderedLevels <<- c(orderNames[!(orderNames %in% orderedLevels)],orderedLevels)
        }
        
        graphData$reorderedCat <- factor(as.character(graphData$category), levels=orderedLevels)
      
      }
        
    graphData
    
  })

# #############################################################    
   
  output$mymap <- leaflet::renderLeaflet({
    
    map <- leaflet(height = "50px") %>%
      addProviderTiles("CartoDB.Positron") %>%
      setView(lng = -83.5, lat = 44.5, zoom=6) 
    
    map
    
  })
  
  observe({
    
    chemGroup <- chemicalSummary()
    
    # chemGroup <- filter(chemGroup, class == "Herbicide")
    
    meanEARlogic <- input$meanEAR
    
    maxEARWords <- ifelse(meanEARlogic,"meanEAR","maxEAR")

    if(input$radioMaxGroup == "1"){
      typeWords <- "groups"
    } else if (input$radioMaxGroup == "2"){
      typeWords <- "chemicals"
    } else {
      typeWords <- "classes"
    }
      
    statsOfGroupOrdered <- statsOfGroupOrdered()
    
    sumStat <- chemGroup %>%
      group_by(site, date) %>%
      summarise(sumEAR = sum(EAR),
                hits=as.numeric(any(hits > 0))) %>%
      data.frame() %>%
      group_by(site) %>%
      summarise(nSamples = n(),
                meanEAR=ifelse(as.logical(meanEARlogic),mean(sumEAR,na.rm=TRUE),max(sumEAR,na.rm=TRUE)),
                freq=sum(hits)/n()) %>%
      data.frame()

    
    siteToFind <- unique(sumStat$site)

    mapData <- right_join(stationINFO[,c("shortName", "Station.Name", "dec.lat.va","dec.long.va")], sumStat, by=c("shortName"="site"))
    mapData <- left_join(mapData, statsOfGroupOrdered, by=c("shortName"="site", "nSamples"="nSamples"))
    
    mapData <- mapData[!is.na(mapData$dec.lat.va),]
    mapData <- unique(mapData)
    
    col_types <- c("darkblue","dodgerblue","green4","gold1","orange","brown","red")

    if(as.logical(meanEARlogic)){
      counts <- mapData$mean
    } else {
      counts <- mapData$max
    }
    
    mapData$meanEAR[is.na(mapData$meanEAR)] <- 0
    counts[is.na(counts)] <- 0
    
    if(nrow(mapData) > 1){
      leg_vals <- unique(as.numeric(quantile(mapData$meanEAR, probs=c(0,0.01,0.1,0.25,0.5,0.75,0.9,.99,1), na.rm=TRUE)))
      pal = colorBin(col_types, mapData$meanEAR, bins = leg_vals)
      rad <-3*seq(1,4,length.out = 16)
      # rad <- 1.5*seq(5000,20000, 1000)
      mapData$sizes <- rad[as.numeric(cut(counts, breaks=16))]
    } else {
      leg_vals <- unique(as.numeric(quantile(c(0,mapData$meanEAR), probs=c(0,0.01,0.1,0.25,0.5,0.75,0.9,.99,1), na.rm=TRUE)))
      pal = colorBin(col_types, c(0,mapData$maxEAR), bins = leg_vals)
      mapData$sizes <- 3
      # mapData$sizes <- 1.5*12000
    }

    # leg_vals <- c(0,0.0001,0.001,0.01,0.1,1,10)
    # pal = colorBin(col_types, c(0,mapData$maxEAR), bins = leg_vals)

    map <- leafletProxy("mymap", data=mapData) %>%
      clearMarkers() %>%
      clearControls() %>%
      addCircleMarkers(lat=~dec.lat.va, lng=~dec.long.va, 
                 popup=paste0('<b>',mapData$Station.Name,"</b><br/><table>",
                              "<tr><td>",maxEARWords,": </td><td>",sprintf("%.1f",mapData$meanEAR),'</td></tr>',
                              "<tr><td>Number of Samples: </td><td>",mapData$nSamples,'</td></tr>',
                              "<tr><td>Frequency: </td><td>",sprintf("%.1f",mapData$freq),'</td></tr>',
                              "<tr><td>Number of ",typeWords," with hits: </td><td>",counts,'</td></tr>',
                              '</table>') ,
                 fillColor = ~pal(meanEAR), 
                 # weight = 1,
                 # color = "black",
                 fillOpacity = 0.8, 
                 radius = ~sizes, 
                 stroke=FALSE,
                 opacity = 0.8) 
    

    if(length(siteToFind) > 1){
      map <- addLegend(map,
        position = 'bottomleft',
        pal=pal,
        values=~meanEAR,
        opacity = 0.8,
        labFormat = labelFormat(digits = 2), #transform = function(x) as.integer(x)),
        title = ifelse(meanEARlogic,'Mean EAR','Max EAR'))
        # title = "Herbicide Max EAR")
    }
    
    map
    
  })
  
############################################################# 
  
  output$TableHeader <- renderUI({
    HTML(paste("<h4>", input$group,"-",input$data,": ",input$sites, "</h4>"))
  })
  
  output$mapFooter <- renderUI({
    
    statsOfGroupOrdered <- statsOfGroupOrdered()
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
  
  output$BoxHeader <- renderUI({
    HTML(paste("<h4>", input$group,"-",input$data,": ",input$sites, "</h4>"))
  })

  output$TableHeaderColumns <- renderUI({
    HTML(paste("<h4>", input$data,": ",input$sites, "</h4>"))
  })
  
  output$siteHitText <- renderUI({
    
    if(input$sites == "All" | input$sites == "Potential 2016"){
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
  
  output$endpointGraph <- renderPlot({ 

    filterBy <- input$epGroup
    meanEARlogic <- as.logical(input$meanEAR)
    hitThres <- hitThresValue()
    
    filterCat <- switch(as.character(input$radioMaxGroup),
                        "1" = "choices",
                        "2" = "chnm",
                        "3" = "class")
    
    boxData <- chemicalSummary()
  
    graphData <- boxData %>%
      # filter(!is.na(category)) %>%
      group_by(site,date,category,endPoint) %>%
      summarise(sumEAR=sum(EAR)) %>%
      data.frame() %>%
      group_by(site, category,endPoint) %>%
      summarise(meanEAR=ifelse(meanEARlogic,mean(sumEAR),max(sumEAR))) %>%
      data.frame() %>%
      mutate(category=as.character(category)) 

    if(filterBy != "All"){
      graphData <- graphData %>%
        filter_(paste0("category == '", filterBy,"'"))
      
      countNonZero <- graphData %>%
        group_by(endPoint) %>%
        summarise(nonZero = as.character(sum(meanEAR>0)),
                  hits = as.character(sum(meanEAR>hitThres)))
      
      countNonZero$hits[countNonZero$hits == "0"] <- ""
      
      namesToPlotEP <<- as.character(countNonZero$endPoint)
      nSamplesEP <<- countNonZero$nonZero
      nHitsEP <<- countNonZero$hits
    }
    
    orderColsBy <- graphData %>%
      group_by(endPoint) %>%
      summarise(median = quantile(meanEAR[meanEAR != 0],0.5)) %>%
      arrange(median)
    
    orderedLevelsEP <<- orderColsBy$endPoint
    
    if(any(is.na(orderColsBy$median))){
      orderedLevelsEP <<- c(orderColsBy$endPoint[is.na(orderColsBy$median)],
                          orderColsBy$endPoint[!is.na(orderColsBy$median)])
    }
    
    graphData$endPoint <- factor(graphData$endPoint, levels = orderedLevelsEP)
    
    stackedPlot <- ggplot(graphData)+
      scale_y_log10(paste(ifelse(meanEARlogic,"Mean","Maximum"), "EAR Per Site"),labels=fancyNumbers) +
      geom_boxplot(aes(x=endPoint, y=meanEAR)) +
      theme_minimal() +
      xlab("") +
      theme(axis.text.y = element_text(vjust = .25,hjust=1)) +
      geom_hline(yintercept = hitThres, linetype="dashed", color="black")
    
    if(filterBy != "All"){
      
      ymin <<- 10^(ggplot_build(stackedPlot)$panel$ranges[[1]]$y.range)[1]
      ymax <<- 10^(ggplot_build(stackedPlot)$panel$ranges[[1]]$y.range)[2]
      
      xmax <<- ggplot_build(stackedPlot)$panel$ranges[[1]]$x.range[2]
      xmin <<- ggplot_build(stackedPlot)$panel$ranges[[1]]$x.range[1]
      
      stackedPlot <- stackedPlot + 
        geom_text(data=data.frame(), aes(x=namesToPlotEP, y=ymin,label=nSamplesEP),size=5)  +
        geom_text(data=data.frame(), aes(x=namesToPlotEP, y=ymax,label=nHitsEP),size=5) 
      
      df1 <- data.frame(y = c(ymin,hitThres,ymax), text = c("# Non Zero","Hit Threshold","# Hits"), stringsAsFactors = FALSE)
      
      for(i in 1:3){
        stackedPlot <- stackedPlot + 
          annotation_custom(
            grob = textGrob(label = df1$text[i], gp = gpar(cex = 0.75)),
            ymin = log10(df1$y[i]),      # Vertical position of the textGrob
            ymax = log10(df1$y[i]),
            xmin = xmax+0.05,  # Note: The grobs are positioned outside the plot area
            xmax = xmax+0.05)
      }
    }
    stackedPlot <- stackedPlot +
      coord_flip()
    
    # Code to override clipping
    stackedPlot <- ggplot_gtable(ggplot_build(stackedPlot))
    stackedPlot$layout$clip[stackedPlot$layout$name == "panel"] <- "off"
    
    print(grid.draw(stackedPlot))
    
  })
  
  output$hitsTable <- DT::renderDataTable({    
    
    boxData <- chemicalSummary()
    meanEARlogic <- as.logical(input$meanEAR)
    hitThres <- hitThresValue()
    
    if(length(unique(boxData$site)) > 1){
      tableData <- boxData %>%
        group_by(site, choices, category, date) %>%
        summarize(sumEAR = sum(EAR)) %>%
        group_by(site, choices, category) %>%
        summarize(meanEAR = ifelse(meanEARlogic, mean(sumEAR),max(sumEAR))) %>%
          # hits = any(hits > 0)) %>% #is a hit when any EAR is greater than 0.1?
        group_by(choices, category) %>%
        summarize(nSites = sum(meanEAR>hitThres)) %>%
        data.frame() 
    } else {
      tableData <- boxData %>%
        group_by(choices, category, date)
      
      tableData <- tableData %>%
        summarise(sumEAR=sum(EAR)) %>%
        data.frame() %>%
        group_by(choices, category) %>%
        summarise(nSites = sum(sumEAR>hitThres))%>%
        data.frame() 
      
    }
    
    if(input$radioMaxGroup != "1"){
      tableData <- tableData %>%
        reshape(idvar="choices",timevar="category", direction="wide") 
      
      names(tableData) <- gsub("nSites.","",names(tableData))
      names(tableData)[1] <- "Groups"
      
      sumOfColumns <- colSums(tableData[-1],na.rm = TRUE)
      orderData <- order(sumOfColumns,decreasing = TRUE) 
      orderData <- orderData[sumOfColumns[orderData] != 0] + 1
      
      tableData <- tableData[,c(1,orderData)]
      colors <- brewer.pal(9,"Blues") #"RdYlBu"
      
      groups <- tableData$Groups
      
      tableData <- tableData[!is.na(groups),-1,drop=FALSE]
      rownames(tableData) <- groups[!is.na(groups)]
      
      cuts <- seq(0,max(as.matrix(tableData),na.rm=TRUE),length.out = 8)
      
      names(tableData)[names(tableData) == "Human Drug, Non Prescription"] <- "Human Drug"
      names(tableData)[names(tableData) == "Flavor/Fragrance"] <- "Flavor / Fragrance"
    } else {
      tableData <- select(tableData, choices, nSites)
      rownames(tableData) <- tableData$choices
      tableData <- tableData[,-1,drop=FALSE]
    }
    
    tableData1 <- DT::datatable(tableData, extensions = 'Buttons',
                                rownames = TRUE,
                                options = list(scrollX = TRUE,
                                               dom = 'Bfrtip',
                                               buttons = 
                                                 list('colvis', list(
                                                   extend = 'collection',
                                                   buttons = list(list(extend='csv',
                                                                       filename = 'siteHits'),
                                                                  list(extend='excel',
                                                                       filename = 'siteHits'),
                                                                  list(extend='pdf',
                                                                       filename= 'siteHits')),
                                                   text = 'Download',
                                                   filename= 'test'
                                                 )),
                                 pageLength = nrow(tableData),
                                 order=list(list(1,'desc'))))
    if(input$radioMaxGroup != "1"){
      for(i in 1:ncol(tableData)){
        tableData1 <- formatStyle(tableData1, columns = names(tableData)[i], 
                    backgroundColor = styleInterval(cuts = cuts,values = colors),
                    color = styleInterval(0.75*max(tableData,na.rm=TRUE),values = c("black","white")),
                    `font-size` = '17px')
      }
    }
    tableData1
    
  })
  
  output$hitsTableEPs <- DT::renderDataTable({

    boxData <- chemicalSummary()
    meanEARlogic <- as.logical(input$meanEAR)
    catType <- as.numeric(input$radioMaxGroup)
    
    hitThres <- hitThresValue()

    fullData_init <- data.frame(Endpoint="",stringsAsFactors = FALSE)
    fullData <- fullData_init
    
    if(length(unique(boxData$site)) > 1){
      
      for(i in unique(boxData$choices)){
        dataSub <- boxData %>%
          filter(choices == i) %>%
          group_by(site, category, endPoint, date) %>%
          summarize(sumEAR = sum(EAR)) %>%
          group_by(site, category, endPoint) %>%
          summarize(meanEAR = ifelse(meanEARlogic, mean(sumEAR),max(sumEAR))) %>%
          group_by(category, endPoint) %>%
          summarize(nSites = sum(meanEAR>hitThres)) %>%
          data.frame() %>%
          arrange(desc(nSites)) %>%
          reshape(idvar="endPoint",timevar="category", direction="wide") 
        
        names(dataSub) <- gsub("nSites.","",names(dataSub))
        names(dataSub)[1] <- "Endpoint"

        if(ncol(dataSub) > 2){
          dataSub <- dataSub[,c(1,1+which(colSums(dataSub[,-1],na.rm = TRUE) != 0))]
        }
        
        if(is.data.frame(dataSub)){
          if(ncol(dataSub) > 2){
            dataSub <- dataSub[(rowSums(dataSub[,-1],na.rm = TRUE) != 0),]
          } else {
            dataSub <- dataSub[which(dataSub[,-1] != 0 ),]
          }
          
          dataSub <- dataSub %>%
            data.frame() %>%
            mutate(Group = i)
          
          fullData <- full_join(fullData,dataSub)

        }
      }
      
    } else {
      
      for(i in unique(boxData$choices)){
        dataSub <- boxData %>%
          filter(choices == i) %>%
          group_by(category, endPoint, date) %>%
          summarise(sumEAR=sum(EAR)) %>%
          data.frame() %>%
          group_by(endPoint, category) %>%
          summarise(nSites = sum(sumEAR>hitThres))%>%
          data.frame() %>%
          arrange(desc(nSites)) %>%
          reshape(idvar="endPoint",timevar="category", direction="wide") 
        
        names(dataSub) <- gsub("nSites.","",names(dataSub))
        names(dataSub)[1] <- "Endpoint"
        
        if(ncol(dataSub) > 2){
          dataSub <- dataSub[,c(1,1+which(colSums(dataSub[,-1],na.rm = TRUE) != 0))]
        }
        
        if(is.data.frame(dataSub)){
          if(ncol(dataSub) > 2){
            dataSub <- dataSub[(rowSums(dataSub[,-1],na.rm = TRUE) != 0),]
          } else {
            dataSub <- dataSub[which(dataSub[,-1] != 0 ),]
          }
          dataSub <- dataSub %>%
            data.frame() %>%
            mutate(Group = i)
          
          fullData <- full_join(fullData,dataSub)
        }
      }
    }
    
    fullData <- fullData[,c("Endpoint","Group",names(fullData)[!(names(fullData) %in% c("Endpoint","Group"))])]

    fullData <- fullData[-1,]
    
    names(fullData) <- gsub("\\."," ",names(fullData))
    names(fullData)[names(fullData) == "Human Drug  Non Prescription"] <- "Human Drug"
    names(fullData)[names(fullData) == "Flavor/Fragrance"] <- "Flavor / Fragrance"
    
    sumOfColumns <- colSums(fullData[c(-1,-2)],na.rm = TRUE)
    orderData <- order(sumOfColumns,decreasing = TRUE) 
    orderData <- orderData[sumOfColumns[orderData] != 0] + 2
    
    fullData <- fullData[,c(1,2,orderData)]
    colors <- brewer.pal(9,"Blues") #"RdYlBu"
    
    groups <- fullData$Groups
    
    if(catType == 2){
      casKey <- select(boxData, chnm, casrn) %>%
        distinct()
      
      hits <- sapply(fullData, function(x) as.character(x))
      
      for(k in 1:nrow(fullData)){
        for(z in 3:ncol(fullData)){
          if(!is.na(fullData[k,z])){
            hits[k,z] <- createLink(ep = fullData$Endpoint[k], 
                                    cas = casKey$casrn[casKey$chnm == names(fullData)[z]], 
                                    hits = fullData[k,z])
          }
        }
      }
      
      fullData <- hits
    }

    fullData <- DT::datatable(fullData, extensions = 'Buttons',
                                escape = FALSE,
                                rownames = FALSE,
                                options = list(dom = 'Bfrtip',
                                               buttons = 
                                                 list('colvis', list(
                                                   extend = 'collection',
                                                   buttons = list(list(extend='csv',
                                                                       filename = 'epHits'),
                                                                  list(extend='excel',
                                                                       filename = 'epHits'),
                                                                  list(extend='pdf',
                                                                       filename= 'epHits')),
                                                   text = 'Download'
                                                 )),
                                               scrollX = TRUE,
                                               pageLength = nrow(fullData),
                                               order=list(list(2,'desc'))))
    
  })

})