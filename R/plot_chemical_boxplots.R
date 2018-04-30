

plot_chemical_boxplots <- function(chemicalSummary, 
                                   manual_remove=NULL,
                                   mean_logic = FALSE,
                                   plot_ND = TRUE,
                                   font_size = NA,
                                   title = NA,
                                   pallette = NA,
                                   hit_threshold = NA){
  
  site <- EAR <- sumEAR <- meanEAR <- groupCol <- nonZero <- ".dplyr"
  chnm <- Class <- maxEAR <- x <- y <- ".dplyr"

  cbValues <- c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628",
                "#DCDA4B","#999999","#00FFFF","#CEA226","#CC79A7","#4E26CE",
                "#FFFF00","#78C15A","#79AEAE","#FF0000","#00FF00","#B1611D",
                "#FFA500","#F4426e", "#800000", "#808000")
  
  if(!plot_ND){
    chemicalSummary <- chemicalSummary[chemicalSummary$EAR > 0,]
  }
  
  if(length(unique(chemicalSummary$Class)) > length(cbValues)){
    n <- length(unique(chemicalSummary$Class))
    
    if(n > 20 & n<30){
      cbValues <- c(brewer.pal(n = 12, name = "Set3"),
                    brewer.pal(n = 8, name = "Set2"),
                    brewer.pal(n = max(c(3,n-20)), name = "Set1"))
    } else if (n <= 20){
      cbValues <- c(brewer.pal(n = 12, name = "Set3"),
                    brewer.pal(n =  max(c(3,n-12)), name = "Set2"))     
    } else {
      cbValues <- colorRampPalette(brewer.pal(11,"Spectral"))(n)
    }

  }
  
  single_site <- length(unique(chemicalSummary$site)) == 1
  
  y_label <- fancyLabels(category = "Chemical", mean_logic, single_site)
  
  if(single_site){
    
    countNonZero <- chemicalSummary %>%
      select(chnm, Class, EAR) %>%
      group_by(chnm, Class) %>%
      summarize(nonZero = as.character(sum(EAR>0)),
                hits = as.character(sum(EAR > hit_threshold)))
    
    countNonZero$hits[countNonZero$hits == "0"] <- ""
    
    label <- "# Endpoints"
    
    pretty_range <- range(chemicalSummary$EAR[chemicalSummary$EAR > 0])
    pretty_logs <- 10^(-10:10)
    log_index <- which(pretty_logs < pretty_range[2] & pretty_logs > pretty_range[1])
    log_index <- c(log_index[1]-1,log_index, log_index[length(log_index)]+1)
    pretty_logs_new <-  pretty_logs[log_index]
    
    toxPlot_All <- ggplot(data=chemicalSummary)
    
    if(!all(is.na(pallette))){
      toxPlot_All <- toxPlot_All +
        geom_boxplot(aes(x=chnm, y=EAR, fill=chnm),lwd=0.1,outlier.size=1) +
        scale_fill_manual(values = pallette) +
        theme(legend.position = "none")
    } else {
      toxPlot_All <- toxPlot_All +
        geom_boxplot(aes(x=chnm, y=EAR, fill=Class),
                     lwd=0.1,outlier.size=1)
    }
    
  } else {
    
    graphData <- graph_chem_data(chemicalSummary=chemicalSummary, 
                                 manual_remove=manual_remove,
                                 mean_logic=mean_logic)
    
    pretty_range <- range(graphData$maxEAR[graphData$maxEAR > 0])
    pretty_logs <- 10^(-10:10)
    log_index <- which(pretty_logs < pretty_range[2] & pretty_logs > pretty_range[1])
    log_index <- c(log_index[1]-1,log_index, log_index[length(log_index)]+1)
    pretty_logs_new <-  pretty_logs[log_index] 
    
    countNonZero <- graphData %>%
      select(chnm, Class, maxEAR) %>%
      group_by(chnm, Class) %>%
      summarize(nonZero = as.character(sum(maxEAR>0)),
                hits = as.character(sum(maxEAR > hit_threshold)))
    
    countNonZero$hits[countNonZero$hits == "0"] <- ""
    
    label <- "# Sites"
    toxPlot_All <- ggplot(data=graphData) 
    
    if(!all(is.na(pallette))){
      toxPlot_All <- toxPlot_All +
        geom_boxplot(aes(x=chnm, y=maxEAR, fill=chnm),lwd=0.1,outlier.size=1) +
        scale_fill_manual(values = pallette) +
        theme(legend.position = "none")
    } else {
      toxPlot_All <- toxPlot_All +
        geom_boxplot(aes(x=chnm, y=maxEAR, fill=Class),
                     lwd=0.1,outlier.size=1)
    }
 
  }
  
  toxPlot_All <- toxPlot_All +
    scale_y_log10(y_label, labels=fancyNumbers,breaks=pretty_logs_new)  +
    theme_bw() +
    scale_x_discrete(drop = TRUE) +
    coord_flip() +
    geom_hline(yintercept = hit_threshold, linetype="dashed", color="black") +
    theme(axis.text = element_text( color = "black"),
          axis.title.y = element_blank(),
          panel.background = element_blank(),
          plot.background = element_rect(fill = "transparent",colour = NA),
          strip.background = element_rect(fill = "transparent",colour = NA),
          strip.text.y = element_blank(),
          panel.border = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(hjust = 0.5))  
  
  if(all(is.na(pallette))){
    toxPlot_All <- toxPlot_All +
      scale_fill_manual(values = cbValues, drop=FALSE) +
      guides(fill=guide_legend(ncol=6)) +
      theme(legend.position="bottom",
            legend.justification = "left",
            legend.background = element_rect(fill = "transparent", colour = "transparent"),
            legend.title=element_blank(),
            legend.key.height = unit(1,"line"))
  }
    
  if(!is.na(font_size)){
    toxPlot_All <- toxPlot_All +
      theme(axis.text = element_text(size = font_size),
            axis.title =   element_text(size=font_size))
  }
  
  plot_info <- ggplot_build(toxPlot_All)
  layout_stuff <- plot_info$layout
  
  if(packageVersion("ggplot2") >= "2.2.1.9000"){
    ymin <- 10^(layout_stuff$panel_scales_y[[1]]$range$range[1])
    ymax <- 10^(layout_stuff$panel_scales_y[[1]]$range$range[2])
  } else {
    ymin <- 10^(layout_stuff$panel_ranges[[1]]$x.range[1])
    ymax <- 10^(layout_stuff$panel_ranges[[1]]$x.range[2])
  }
  
  toxPlot_All_withLabels <- toxPlot_All +
    geom_text(data=countNonZero, aes(x=chnm,label=nonZero, y=ymin), size = ifelse(is.na(font_size),2,0.30*font_size)) +
    geom_text(data=data.frame(x = Inf, y=ymin, label = label, stringsAsFactors = FALSE), 
            aes(x=x,  y=y, label = label),
            size=ifelse(is.na(font_size),3,0.30*font_size)) 
  
  nHitsEP <- countNonZero$hits
  
  if(isTRUE(sum(as.numeric(nHitsEP), na.rm = TRUE) > 0)) {
    toxPlot_All_withLabels <- toxPlot_All_withLabels +
      geom_text(data=countNonZero, aes(x=chnm, y=ymax,label=nHitsEP),size=ifelse(is.na(font_size),3,0.30*font_size)) +
      geom_text(data=data.frame(x = Inf, y=ymax, label = "# Hits", stringsAsFactors = FALSE), 
                aes(x = x,  y=y, label = label),
                size=ifelse(is.na(font_size),3,0.30*font_size))
  }
  
  if(!all(is.na(title))){
    toxPlot_All_withLabels <- toxPlot_All_withLabels +
      ggtitle(title)
    
    if(!is.na(font_size)){
      toxPlot_All_withLabels <- toxPlot_All_withLabels +
        theme(plot.title = element_text(size=font_size))
    }
  }
  
  if(!is.na(hit_threshold)) {
    toxPlot_All_withLabels <- toxPlot_All_withLabels +
      geom_text(data=data.frame(x = Inf, y=hit_threshold, label = "Threshold", stringsAsFactors = FALSE), 
                aes(x = x,  y=y, label = label),
                size=ifelse(is.na(font_size),3,0.30*font_size))
  }
  
  return(toxPlot_All_withLabels)
  
}

#' @export
#' @rdname plot_tox_boxplots
graph_chem_data <- function(chemicalSummary, 
                            manual_remove=NULL,
                            mean_logic = FALSE){
  
  site <- chnm <- Class <- EAR <- sumEAR <- maxEAR <- ".dplyr"

  mean_logic <- as.character(mean_logic)
  match.arg(mean_logic, c("mean","max","noSum","TRUE","FALSE"))
  
  if(mean_logic %in% c("TRUE","mean")){
    mean_logic <- TRUE
  }
  if(mean_logic %in% c("FALSE","max")){
    mean_logic <- FALSE
  }
  
  if(mean_logic == "noSum"){
    graphData <- chemicalSummary %>%
      group_by(site, chnm, Class) %>%
      summarise(maxEAR=max(EAR)) %>%
      data.frame()     
  } else {
    graphData <- chemicalSummary %>%
      group_by(site,date,chnm, Class) %>%
      summarise(sumEAR=sum(EAR)) %>%
      data.frame() %>%
      group_by(site, chnm, Class) %>%
      summarise(maxEAR=ifelse(mean_logic,mean(sumEAR),max(sumEAR))) %>%
      data.frame() 
  }
  
  if(!is.null(manual_remove)){
    graphData <- filter(graphData, !(chnm %in% manual_remove))
  }
  
  return(graphData)
}

