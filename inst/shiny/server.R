library(dplyr)
library(ggplot2)
library(DT)
library(leaflet)
library(data.table)
library(tidyr)
library(toxEval)
library(RColorBrewer)
library(grid)

endPointInfo <- endPointInfo

endPointInfo <- endPointInfo[,c(39,51,30)]
endPointInfo <- endPointInfo[-grep("BSK",endPointInfo$assay_component_endpoint_name),]
endPointInfo <- endPointInfo[-grep("APR",endPointInfo$assay_component_endpoint_name),]

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

interl3 <- function (a,b,cv) {
  n <- min(length(a),length(b),length(cv))
  p1 <- as.vector(rbind(a[1:n],b[1:n],cv[1:n]))
  p2 <- c(a[-(1:n)],b[-(1:n)],cv[-(1:n)])
  c(p1,p2)
}

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

makePlots <- function(boxData, noLegend, boxPlot){
  
  siteToFind <- unique(boxData$site)

  if(length(siteToFind) > 1){
    graphData <- boxData %>%
      filter(!is.na(category)) %>%
      group_by(site,date,category) %>%
      summarise(sumEAR=sum(EAR)) %>%
      data.frame() %>%
      group_by(site, category) %>%
      summarise(meanEAR=mean(sumEAR)) %>%
      data.frame() %>%
      mutate(category=as.character(category)) 
    
    orderColsBy <- graphData %>%
      group_by(category) %>%
      summarise(median = quantile(meanEAR[meanEAR != 0],0.5)) %>%
      arrange(median)
    
    orderedLevels <<- orderColsBy$category[!is.na(orderColsBy$median)] #The !is.na takes out any category that was all censo

    graphData$reorderedCat <- factor(as.character(graphData$category), levels=orderedLevels)
    graphData <- filter(graphData, !is.na(reorderedCat))
    
    lowerPlot <- ggplot(graphData)+
      scale_y_log10("Mean EAR Per Site",labels=fancyNumbers) # fance labels
    labelsText <<- "nSites"
    
    countNonZero <- graphData %>%
      group_by(category) %>%
      summarise(nonZero = as.character(sum(meanEAR>0)),
                hits = as.character(sum(meanEAR>0.1)))
    countNonZero$nonZero[countNonZero$nonZero == "0"] <- ""
    countNonZero$hits[countNonZero$hits == "0"] <- ""
    
    namesToPlot <<- as.character(countNonZero$category)
    nSamples <<- countNonZero$nonZero
    nHits <<- countNonZero$hits
    
  } else {
    
    graphData <- boxData %>%
      group_by(site,date,category) %>%
      summarise(sumEAR=sum(EAR)) %>%
      data.frame() %>%
      rename(meanEAR=sumEAR)%>%
      mutate(category=factor(category)) %>%
      filter(!is.na(category)) 
    
    orderColsBy <- graphData %>%
      mutate(category = as.character(category)) %>%
      group_by(category) %>%
      summarise(median = median(meanEAR[meanEAR != 0])) %>%
      arrange(!is.na(median),median)
    
    orderedLevels <<- orderColsBy$category#[!is.na(orderColsBy$median)] #The !is.na takes out any category that was all censo
    graphData$reorderedCat <- factor(as.character(graphData$category), levels=orderedLevels)
    graphData <- filter(graphData, !is.na(reorderedCat))
    
    lowerPlot <- ggplot(graphData)+
      scale_y_log10("Sum of EAR",labels=fancyNumbers) # labels=fancyNumbers

    labelsText <<- "nSamples"
    
    countNonZero <- graphData %>%
      group_by(category) %>%
      summarise(nonZero = as.character(sum(meanEAR>0)),
                hits = as.character(sum(meanEAR>0.1)))
    countNonZero$nonZero[countNonZero$nonZero == "0"] <- ""
    countNonZero$hits[countNonZero$hits == "0"] <- ""
    
    namesToPlot <<- as.character(countNonZero$category)
    nSamples <<- countNonZero$nonZero
    nHits <<- countNonZero$hits

  }
  
  if(!boxPlot){
    lowerPlot <- lowerPlot + geom_point(aes(x=reorderedCat, y=meanEAR, color=category, size=3))
  } else {
    lowerPlot <- lowerPlot + geom_boxplot(aes(x=reorderedCat, y=meanEAR, fill=category)) 
  }
  
  text.size <<- ifelse(length(orderedLevels) > 15, 13, 16)
  
  lowerPlot <- lowerPlot + 
    theme(legend.position = "none") +
    theme(axis.text.y = element_text(size=text.size,color = "black"), 
          axis.text.x = element_text(size=16,color = "black", vjust = 0),
          axis.title = element_text(size=16)) +
    xlab("") +
    geom_hline(yintercept = 0.1, linetype="dashed", color="red")  +
    scale_x_discrete(drop=FALSE) +
    scale_fill_discrete(drop=FALSE)

  ymin <<- 10^(ggplot_build(lowerPlot)$panel$ranges[[1]]$y.range)[1]
  ymax <<- 10^(ggplot_build(lowerPlot)$panel$ranges[[1]]$y.range)[2]
  
  xmax <<- ggplot_build(lowerPlot)$panel$ranges[[1]]$x.range[2]
  xmin <<- ggplot_build(lowerPlot)$panel$ranges[[1]]$x.range[1]

  lowerPlot <- lowerPlot + 
    geom_text(data=data.frame(), aes(x=namesToPlot, y=ymin,label=nSamples),size=5)  +
    geom_text(data=data.frame(), aes(x=namesToPlot, y=ymax,label=nHits),size=5) 
  
  df1 <- data.frame(y = c(ymin,0.1,ymax), text = c("# Non Zero","Hit Threshold","# Hits"), stringsAsFactors = FALSE)
  
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
    # scale_x_continuous(limits = c(xmin,xmax+0.5)) + 
    coord_flip() 
  
  # Code to override clipping
  lowerPlot <- ggplot_gtable(ggplot_build(lowerPlot))
  lowerPlot$layout$clip[lowerPlot$layout$name == "panel"] <- "off"

  siteLimits <- stationINFO %>%
    filter(shortName %in% unique(boxData$site))
  
  if(length(siteToFind) > 1){

    graphData <- graphData %>%
      arrange(as.character(category)) 
    
    if(all(siteLimits$shortName %in% sitesOrdered)){
      siteLimits <- mutate(siteLimits, shortName = factor(shortName, levels=sitesOrdered[sitesOrdered %in% siteLimits$shortName]))
    } else {
      siteLimits <- mutate(siteLimits, shortName = factor(shortName))
    }
    
    upperPlot <- ggplot(graphData, aes(x=site, y=meanEAR, fill = category)) +
      geom_bar(stat="identity") +
      # facet_wrap(~stat, nrow=2, ncol=1, scales = "free_y") + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.25, 
                                       colour=siteLimits$lakeColor)) +
      scale_x_discrete(limits=levels(siteLimits$shortName),drop=FALSE) +
      xlab("") +
      ylab("Mean EAR per Site") +
      scale_fill_discrete("", drop=FALSE) 
    
    if(noLegend){
      upperPlot <- upperPlot + 
      guides(fill=FALSE) 
    }

  } else {
    
    chemGroupBPOneSite <- boxData %>%
      select(-site)  %>%
      arrange(as.character(category)) 
    
    upperPlot <- ggplot(chemGroupBPOneSite, aes(x=date, y=EAR, fill = category)) +
      geom_bar(stat="identity")+
      theme(axis.text.x=element_blank(),
            axis.ticks=element_blank())+
      xlab("Individual Samples") + 
      ylab("EAR") +
      scale_fill_discrete("", drop=FALSE) +
      labs(fill="") 
      # scale_y_log10("EAR")
    
    if(noLegend){
      upperPlot <- upperPlot +
      guides(fill=FALSE) 
    } 
  }
  
  return(list(upperPlot=upperPlot, lowerPlot=lowerPlot))
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

  observe({
    
    columnName <- input$groupCol
    
    updateSelectInput(session, "group", selected = "All")
    
  })
  
  observe({
       
       columnName <- input$groupCol
       orderBy <- endPointInfo[,columnName]
       orderNames <- names(table(orderBy))
       nEndPoints <- as.integer(table(orderBy))
       df <- data.frame(orderNames,nEndPoints,stringsAsFactors = FALSE) %>%
         arrange(desc(nEndPoints))
       
       if(length(grep("background",df$orderNames)) > 0){
         df <- df[-grep("background",df$orderNames),]
       }
       
       if(length(grep("cell cycle",df$orderNames)) > 0){
         df <- df[-grep("cell cycle",df$orderNames),]
       }
       
       dropDownHeader <- c("All",paste0(df$orderNames," (",df$nEndPoints,")"))
        
       updateSelectInput(session, "group",
        choices = setNames(c("All",df$orderNames),dropDownHeader)
       )
    })
  

      
#############################################################   
    chemicalSummary <- reactive({
      
      path <- pathToApp
      groupCol <- input$groupCol
      group <- input$group
      radioMaxGroup <- input$radioMaxGroup
      
      if (input$data == "Water Sample"){
        chemicalSummary <- readRDS(file.path(path,"chemicalSummaryV2.rds"))
        stationINFO <<- readRDS(file.path(path,"sitesOWC.rds"))
      } else if (input$data == "Passive Samples"){
        chemicalSummary <- readRDS(file.path(path,"chemicalSummaryPassive.rds"))
        chemicalSummary$date <- rep(as.POSIXct(as.Date("1970-01-01")),nrow(chemicalSummary))
        stationINFO <<- readRDS(file.path(path,"sitesOWC.rds"))
      } else if (input$data == "Duluth"){
        chemicalSummary <- readRDS(file.path(path,"chemSummeryDL.rds"))
        stationINFO <<- readRDS(file.path(path,"sitesDuluth.rds"))
      } else if (input$data == "NPS"){
        chemicalSummary <- readRDS(file.path(path,"chemNPS.rds"))
        stationINFO <<- readRDS(file.path(path,"npsSite.rds"))
      }
      
      
      ep <- data.frame(endPointInfo[,c("assay_component_endpoint_name", groupCol)])
      if(length(grep("background",ep[,2])) > 0){
        ep <- ep[-grep("background",ep[,2]),]
      }
      if(length(grep("cell cycle",ep[,2])) > 0){
        ep <- ep[-grep("cell cycle",ep[,2]),]
      }
      ep <- ep[!is.na(ep[,2]),]
      

      
      chemicalSummary <- rename(chemicalSummary, assay_component_endpoint_name=endPoint)%>%
        filter(assay_component_endpoint_name %in% ep$assay_component_endpoint_name ) %>%
        data.table() %>%
        left_join(data.table(ep)) %>%
        data.frame() %>%
        rename(endPoint = assay_component_endpoint_name) %>%
        select_("hits","EAR","chnm","class","date","choices"=groupCol,"site","endPoint","endPointValue") %>%
        left_join(stationINFO[,c("fullSiteID","shortName")], by=c("site"="fullSiteID")) %>%
        select(-site) %>%
        rename(site=shortName)
      
      if(group != "All"){
        endPointInfoSub <- select_(endPointInfo, "assay_component_endpoint_name", groupCol) %>%
          distinct() %>%
          rename(endPoint=assay_component_endpoint_name) %>%
          filter_(paste0(groupCol," =='", group, "'"))
        
        chemicalSummary <- filter(chemicalSummary, endPoint %in% endPointInfoSub$endPoint)
      }
    
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
        summarise(maxEAR = max(sumEAR),
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

      siteToFind <- unique(statsOfGroup$site)
      
      statsOfGroupOrdered <- statsOfGroup %>%
        group_by(site, date,category) %>%
        summarise(sumEAR = sum(EAR)) %>%
        group_by(site,category) %>%
        summarise(max = sum(sumEAR > 0.1),
                  mean = sum(mean(sumEAR) > 0.1),
                  hit = as.numeric(any(sumEAR > 0.1)),
                  nSamples = n())%>%
        data.frame()

      if(length(siteToFind) > 1){
        statsOfGroupOrdered <- statsOfGroupOrdered %>%
          group_by(site) %>%
          summarise(max=sum(max > 0),
                    mean=sum(mean > 0),
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

        tableGroup <- DT::datatable(statsOfGroupOrdered, 
                                      rownames = FALSE,
                                      options = list(dom = 'ft',
                                                     pageLength = nrow(statsOfGroupOrdered), 
                                                     order=list(list(colToSort,'desc'))))

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
          
          tableGroup <- DT::datatable(statsOfGroupOrdered[,c("category","max","nSamples")], 
                                      rownames = FALSE,
                                      colnames = c('hits' = 2),
                                      options = list(dom = 'ft',
                                                     pageLength = nrow(statsOfGroupOrdered), 
                                                     order=list(list(1,'desc'))))

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
      
      colors <- brewer.pal(length(maxEARS),"Blues") #"RdYlBu"
      tableSumm <- DT::datatable(statCol, 
                                 rownames = FALSE,
                                 options = list(dom = 'ft',
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
      print(groupPlots()$upperPlot)
    })

    output$graphGroup <- renderPlot({ 
      print(grid.draw(groupPlots()$lowerPlot))
    })

    groupPlots <- reactive({
      
      boxData <- chemicalSummary()
      
      noLegend <- TRUE
      boxGraph <- TRUE
      siteToFind <- unique(boxData$site)

      radioMaxGroup <- input$radioMaxGroup
      boxGraph <-  !(input$data == "Passive Samples" & length(siteToFind) == 1) 

      return(makePlots(boxData, noLegend, boxGraph))
    })

# #############################################################    
#     
      output$mymap <- leaflet::renderLeaflet({
        
        leaflet(height = "50px") %>%
          # addProviderTiles("Esri.WorldPhysical") %>%
          addProviderTiles("CartoDB.Positron") %>%
          setView(lng = -83.5, lat = 44.5, zoom=6) 
        
      })
      
      observe({
        
        chemGroup <- chemicalSummary()
        
        if (input$data == "Passive Samples"){
          sumStat <- chemGroup %>%
            group_by(site) %>%
            summarise(meanEAR = mean(EAR)) %>%
            mutate(nSamples = 1) %>%
            mutate(freq = NA)        
        } else {
          sumStat <- chemGroup %>%
            group_by(site, date) %>%
            summarise(sumEAR = sum(EAR),
                      hits=as.numeric(any(hits > 0))) %>%
            data.frame() %>%
            group_by(site) %>%
            summarise(nSamples = n(),
                      meanEAR=mean(sumEAR,na.rm=TRUE),
                      freq=sum(hits)/n()) %>%
            data.frame()
        } 
        
        siteToFind <- unique(sumStat$site)

        mapData <- right_join(stationINFO[,c("shortName", "Station.Name", "dec.lat.va","dec.long.va")], sumStat, by=c("shortName"="site"))
        mapData <- mapData[!is.na(mapData$dec.lat.va),]
        
        col_types <- c("darkblue","dodgerblue","green4","gold1","orange","brown","red")
  
        if(nrow(mapData) > 1){
          leg_vals <- unique(as.numeric(quantile(mapData$meanEAR, probs=c(0,0.01,0.1,0.25,0.5,0.75,0.9,.99,1), na.rm=TRUE)))
          pal = colorBin(col_types, mapData$maxEAR, bins = leg_vals)
          rad <-3*seq(1,4,length.out = 16)
          # rad <- 1.5*seq(5000,20000, 1000)
          mapData$sizes <- rad[as.numeric(cut(mapData$nSamples, breaks=16))]
        } else {
          leg_vals <- unique(as.numeric(quantile(c(0,mapData$meanEAR), probs=c(0,0.01,0.1,0.25,0.5,0.75,0.9,.99,1), na.rm=TRUE)))
          pal = colorBin(col_types, c(0,mapData$maxEAR), bins = leg_vals)
          mapData$sizes <- 3
          # mapData$sizes <- 1.5*12000
        }
  
        # pal = colorBin(col_types, c(0,355), bins = c(0,0.1,1,7,60,250,355))
        
        map <- leafletProxy("mymap", data=mapData) %>%
          clearMarkers() %>%
          clearControls() %>%
          addCircleMarkers(lat=~dec.lat.va, lng=~dec.long.va, 
                     popup=paste0('<b>',mapData$Station.Name,"</b><br/><table>",
                                  "<tr><td>maxEAR: </td><td>",sprintf("%.1f",mapData$maxEAR),'</td></tr>',
                                  "<tr><td>Number of Samples: </td><td>",mapData$nSamples,'</td></tr>',
                                  "<tr><td>Frequency: </td><td>",sprintf("%.1f",mapData$freq),'</td></tr></table>') ,
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
            labFormat = labelFormat(digits = 1), #transform = function(x) as.integer(x)),
            title = 'Mean EAR')
        }
        
        map
        
      })
      
  ############################################################# 
      
      output$TableHeader <- renderUI({
        HTML(paste("<h4>", input$group,"-",input$data,": ",input$sites, "</h4>"))
      })
      
      output$mapFooter <- renderUI({
        if(input$data == "Water Sample"){
          HTML("<h5>Size range represents number of collected samples from 1-64</h5>")
        } else if (input$data == "Duluth") {
          HTML("<h5>Size range represents number of collected samples from 1-26</h5>")
        } else {
          HTML("<h5>One sample per site</h5>")
        }
        
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
        
        if(input$sites == "All"){
          siteToFind <- summaryFile$site
        } else if (input$sites == "Potential 2016"){
          siteToFind <- df2016$Site
        } else {
          siteToFind <- input$sites
        }
        
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
        
        filterCat <- switch(as.character(input$radioMaxGroup),
                            "1" = "choices",
                            "2" = "chnm",
                            "3" = "class")
        
        boxData <- chemicalSummary()
      
        boxData3 <- boxData %>%
          filter(endPoint %in% unique(endPoint[EAR > 0.1]))
  
        if(filterBy != "All"){
          boxData3 <- boxData3 %>%
            filter_(paste0(filterCat," == '", filterBy,"'"))
        }
        
        stackedPlot <- ggplot(boxData3, aes(x = cat, y = EAR))+
          scale_y_log10("Mean EAR Per Site") +
          geom_boxplot(aes(x=endPoint, y=EAR)) +
          coord_flip() +
          xlab("") +
          theme(axis.text.y = element_text(vjust = .25,hjust=1)) +
          geom_hline(yintercept = 0.1, linetype="dashed", color="red")
        
        print(stackedPlot)
        
      })
      
      output$hitsTable <- DT::renderDataTable({    
        
        boxData <- chemicalSummary()
        
        if(length(unique(boxData$site)) > 1){
          tableData <- boxData %>%
            group_by(site, choices, category, date) %>%
            summarize(sumEAR = sum(EAR)) %>%
            group_by(site, choices, category) %>%
            summarize(meanEAR = mean(sumEAR)) %>%
              # hits = any(hits > 0)) %>% #is a hit when any EAR is greater than 0.1?
            group_by(choices, category) %>%
            summarize(nSites = sum(meanEAR>0.1)) %>%
            data.frame() 
        } else {
          tableData <- boxData %>%
            group_by(choices, category, date)
          
          tableData <- tableData %>%
            summarise(sumEAR=sum(EAR)) %>%
            data.frame() %>%
            group_by(choices, category) %>%
            summarise(nSites = sum(sumEAR>0.1))%>% #or is a hit when the sum is greater than 0.1?
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
          
          tableData <- tableData[!is.na(groups),-1]
          rownames(tableData) <- groups[!is.na(groups)]
          
          cuts <- seq(0,max(as.matrix(tableData),na.rm=TRUE),length.out = 8)
          
          names(tableData)[names(tableData) == "Human Drug, Non Prescription"] <- "Human Drug"
          names(tableData)[names(tableData) == "Flavor/Fragrance"] <- "Flavor / Fragrance"
        } else {
          tableData <- select(tableData, choices, nSites)
          rownames(tableData) <- tableData$choices
          tableData <- tableData[,-1,drop=FALSE]
        }
        
        tableData1 <- DT::datatable(tableData, # extensions = 'TableTools',
                                    rownames = TRUE,
                                    options = list(dom = 't',
                                                   scrollX = TRUE,
                                                   # initComplete = JS(
                                                   #   "function(settings, json) {",
                                                   #   "$(this.api().table().header()).css({'font-size': 'large'});",
                                                   #   "}"),
                                                   pageLength = nrow(tableData),
                                                   # tableTools = list(sSwfPath = copySWF()),
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

})