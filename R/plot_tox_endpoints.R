#' EndPoint boxplots
#' 
#' The \code{plot_tox_endpoints} function creates a set of boxplots representing EAR 
#' values for each endPoint based on the selected data. A subset of data is first
#' chosen by specifying a group in the filterBy argument. The 
#' filterBy argument must match one of the unique options in the category. 
#' For example, if the category is "Chemical Class", then the filterBy argument 
#' must be one of the defined "Chemical Class" options such as "Herbicide".

#' A boxplot is generated for each endPoint. The EAR values that are used to 
#' create the boxplots are the mean or maximum (as defined by mean_logic) for each 
#' site as described in "Summarizing the data"in the Introduction vignette: 
#' \href{../doc/Introduction.html#summarize_data}{\code{vignette("Introduction", package = "toxEval")}}.
#' 
#' Box plots are standard Tukey representations. See "Box plot details" in the Basic Workflow vignette: 
#' \href{../doc/basicWorkflow.html#box_plot_details}{\code{vignette("basicWorkflow", package = "toxEval")}}
#' for more information.
#' 
#' @param chemical_summary Data frame from \code{\link{get_chemical_summary}}.
#' @param category Either "Biological", "Chemical Class", or "Chemical".
#' @param filterBy Character. Either "All" or one of the filtered categories.
#' @param manual_remove Vector of categories to remove.
#' @param mean_logic Logical.  \code{TRUE} displays the mean sample from each site,
#' \code{FALSE} displays the maximum sample from each site.
#' @param sum_logic logical. \code{TRUE} sums the EARs in a specified grouping,
#' \code{FALSE} does not. \code{FALSE} may be better for traditional benchmarks as
#' opposed to ToxCast benchmarks.
#' @param hit_threshold Numeric threshold defining a "hit".
#' @param font_size Numeric to adjust the axis font size.
#' @param title Character title for plot. 
#' @param palette Vector of color palette for fill. Can be a named vector
#' to specify specific color for specific categories. 
#' @param top_num Integer number of endpoints to include in the graph. If NA, all 
#' endpoints will be included.
#' @export
#' @import ggplot2
#' @importFrom stats median
#' @examples
#' # This is the example workflow:
#' path_to_tox <-  system.file("extdata", package="toxEval")
#' file_name <- "OWC_data_fromSup.xlsx"
#'
#' full_path <- file.path(path_to_tox, file_name)
#' 
#' 
#' tox_list <- create_toxEval(full_path)
#' ACC <- get_ACC(tox_list$chem_info$CAS)
#' ACC <- remove_flags(ACC)
#' 
#' cleaned_ep <- clean_endPoint_info(end_point_info)
#' filtered_ep <- filter_groups(cleaned_ep)
#' chemical_summary <- get_chemical_summary(tox_list, ACC, filtered_ep)
#'
#' plot_tox_endpoints(chemical_summary, filterBy = "Cell Cycle")
#' plot_tox_endpoints(chemical_summary, category = "Chemical Class", filterBy = "PAHs")
#' plot_tox_endpoints(chemical_summary, category = "Chemical", filterBy = "Atrazine")
#' plot_tox_endpoints(chemical_summary, category = "Chemical", 
#'                    filterBy = "Atrazine", top_num = 10)
plot_tox_endpoints <- function(chemical_summary, 
                              category = "Biological",
                              filterBy = "All",
                              manual_remove = NULL,
                              hit_threshold = NA,
                              mean_logic = FALSE, 
                              sum_logic = TRUE,
                              font_size = NA,
                              title = NA,
                              palette = NA, 
                              top_num = NA){
  
  match.arg(category, c("Biological","Chemical Class","Chemical"))

  site <- endPoint <- EAR <- sumEAR <- meanEAR <- x <- y <- ".dplyr"
  
  if(category == "Biological"){
    chemical_summary$category <- chemical_summary$Bio_category
  } else if(category == "Chemical Class") {
    chemical_summary$category <- chemical_summary$Class
  } else {
    chemical_summary$category <- chemical_summary$chnm
  }
      
  single_site <- length(unique(chemical_summary$site)) == 1
  
  if(filterBy != "All"){
    if(!(filterBy %in% unique(chemical_summary$category))){
      stop("filterBy argument doesn't match data")
    }
    chemical_summary <- chemical_summary[chemical_summary["category"] == filterBy,]
  }
  
  y_label <- fancyLabels(category, mean_logic, sum_logic, single_site)
  
  if(single_site){
    
    countNonZero <- chemical_summary %>%
      dplyr::group_by(endPoint) %>%
      dplyr::summarise(nonZero = as.character(sum(EAR>0)),
                hits = as.character(sum(EAR > hit_threshold)))

    countNonZero$hits[countNonZero$hits == "0"] <- ""

    namesToPlotEP <- as.character(countNonZero$endPoint)
    
    nSamplesEP <- countNonZero$nonZero
    nHitsEP <- countNonZero$hits
    
    orderColsBy <- chemical_summary %>%
      dplyr::group_by(endPoint) %>%
      dplyr::summarise(median = mean(EAR[EAR != 0], na.rm = TRUE)) %>%
      dplyr::arrange(median)
    
    orderedLevelsEP <- orderColsBy$endPoint
    
    if(any(is.na(orderColsBy$median))){
      orderedLevelsEP <- c(orderColsBy$endPoint[is.na(orderColsBy$median)],
                           orderColsBy$endPoint[!is.na(orderColsBy$median)])
    }
    
    if(!is.na(top_num)){
      orderedLevelsEP <- orderedLevelsEP[(length(orderedLevelsEP)-top_num+1):length(orderedLevelsEP)]
      chemical_summary <- chemical_summary[chemical_summary[["endPoint"]] %in% orderedLevelsEP,]
      nSamplesEP <- nSamplesEP[(length(nSamplesEP)-top_num+1):length(nSamplesEP)]
      namesToPlotEP <- orderedLevelsEP
    }
    
    pretty_logs_new <-  prettyLogs(chemical_summary$EAR)
    
    chemical_summary$endPoint <- factor(chemical_summary$endPoint, 
                                        levels = orderedLevelsEP)
    
    stackedPlot <- ggplot(data = chemical_summary)+
      scale_y_log10(y_label,labels=fancyNumbers,breaks=pretty_logs_new) +
      theme_minimal() +
      xlab("") +
      theme(axis.text.y = element_text(vjust = .25,hjust=1)) +
      geom_hline(yintercept = hit_threshold, linetype="dashed", color="black")
    
    if(!all(is.na(palette))){
      stackedPlot <- stackedPlot +
        geom_boxplot(aes(x=endPoint, y=EAR, fill = endPoint)) +
        scale_fill_manual(values = palette) +
        theme(legend.position = "none")
    } else {
      stackedPlot <- stackedPlot +
        geom_boxplot(aes(x=endPoint, y=EAR), fill = "steelblue") 
    }
    
  } else {
    
    if(!sum_logic){
      graphData <- chemical_summary %>%
        dplyr::group_by(site, category,endPoint) %>%
        dplyr::summarise(meanEAR=ifelse(mean_logic,mean(EAR),max(EAR))) %>%
        data.frame() %>%
        dplyr::mutate(category=as.character(category))      
    } else {
      
      graphData <- chemical_summary %>%
        dplyr::group_by(site,date,category,endPoint) %>%
        dplyr::summarise(sumEAR=sum(EAR)) %>%
        data.frame() %>%
        dplyr::group_by(site, category,endPoint) %>%
        dplyr::summarise(meanEAR=ifelse(mean_logic,mean(sumEAR),max(sumEAR))) %>%
        data.frame() %>%
        dplyr::mutate(category=as.character(category))      
    }

    pretty_logs_new <- prettyLogs(graphData$meanEAR)
    
    countNonZero <- graphData %>%
      dplyr::group_by(endPoint) %>%
      dplyr::summarise(nonZero = as.character(sum(meanEAR>0)),
                hits = as.character(sum(meanEAR > hit_threshold)))

    countNonZero$hits[countNonZero$hits == "0"] <- ""

    namesToPlotEP <- as.character(countNonZero$endPoint)
    nSamplesEP <- countNonZero$nonZero
    nHitsEP <- countNonZero$hits

    orderColsBy <- graphData %>%
      dplyr::group_by(endPoint) %>%
      dplyr::summarise(median = quantile(meanEAR[meanEAR != 0],0.5)) %>%
      dplyr::arrange(median)
  
    orderedLevelsEP <- orderColsBy$endPoint
  
    if(any(is.na(orderColsBy$median))){
      orderedLevelsEP <- c(orderColsBy$endPoint[is.na(orderColsBy$median)],
                          orderColsBy$endPoint[!is.na(orderColsBy$median)])
    }
  
    if(!is.na(top_num)){
      orderedLevelsEP <- orderedLevelsEP[(length(orderedLevelsEP)-top_num+1):length(orderedLevelsEP)]
      graphData <- graphData[graphData[["endPoint"]] %in% orderedLevelsEP,]
      nSamplesEP <- nSamplesEP[(length(nSamplesEP)-top_num+1):length(nSamplesEP)]
      namesToPlotEP <- orderedLevelsEP
    }
    
    graphData$endPoint <- factor(graphData$endPoint, levels = orderedLevelsEP)
    
    stackedPlot <- ggplot(graphData)+
      scale_y_log10(y_label,labels=fancyNumbers,breaks=pretty_logs_new) +
      theme_minimal() +
      xlab("") +
      theme(axis.text.y = element_text(vjust = .25,hjust=1)) 
    
    if(!is.na(hit_threshold)){
      stackedPlot <- stackedPlot +
        geom_hline(yintercept = hit_threshold, linetype="dashed", color="black")
    }
    
    if(!all(is.na(palette))){
      stackedPlot <- stackedPlot +
        geom_boxplot(aes(x=endPoint, y=meanEAR, fill = endPoint)) +
        scale_fill_manual(values = palette) +
        theme(legend.position = "none")
    } else {
      stackedPlot <- stackedPlot +
        geom_boxplot(aes(x=endPoint, y=meanEAR), fill = "steelblue") 
    }
    
  }

  plot_layout <- ggplot_build(stackedPlot)$layout    
  
  if(utils::packageVersion("ggplot2") >= "2.2.1.9000"){
    ymin <- 10^(plot_layout$panel_scales_y[[1]]$range$range[1])
    ymax <- 10^(plot_layout$panel_scales_y[[1]]$range$range[2])
    
    xmax <- plot_layout$plot_layout$panel_scales_x[[1]]$range$range[1]
    xmin <- plot_layout$plot_layout$panel_scales_x[[1]]$range$range[2]
  } else {
    ymin <- 10^(plot_layout$panel_ranges[[1]]$y.range)[1]
    ymax <- 10^(plot_layout$panel_ranges[[1]]$y.range)[2]

    xmax <- plot_layout$panel_ranges[[1]]$x.range[2]
    xmin <- plot_layout$panel_ranges[[1]]$x.range[1]
  }
  
  label <- "# Sites"
  
  if(single_site){
    label <- "# Chemicals"
  }
  
  stackedPlot <- stackedPlot +
    geom_text(data=data.frame(), aes(x=namesToPlotEP, y=ymin,label=nSamplesEP),size=ifelse(is.na(font_size),3,0.30*font_size)) +
    geom_text(data=data.frame(x = Inf, y=ymin, label = label, stringsAsFactors = FALSE), 
              aes(x = x,  y=y, label = label),
              size=ifelse(is.na(font_size),3,0.30*font_size))     
  
  if(isTRUE(sum(as.numeric(nHitsEP), na.rm = TRUE) > 0)) {
    stackedPlot <- stackedPlot +
      geom_text(data=data.frame(), aes(x=namesToPlotEP, y=ymax,label=nHitsEP),size=ifelse(is.na(font_size),3,0.30*font_size)) +
      geom_text(data=data.frame(x = Inf, y=ymax, label = "# Hits", stringsAsFactors = FALSE), 
              aes(x = x,  y=y, label = label),
              size=ifelse(is.na(font_size),3,0.30*font_size))
  }
  
  if(!is.na(hit_threshold)) {
    stackedPlot <- stackedPlot +
      geom_text(data=data.frame(x = Inf, y=hit_threshold, label = "Threshold", stringsAsFactors = FALSE), 
                aes(x = x,  y=y, label = label),
                size=ifelse(is.na(font_size),3,0.30*font_size))
  }

  if(!is.na(font_size)){
    stackedPlot <- stackedPlot +
      theme(axis.text = element_text(size = font_size),
            axis.title =   element_text(size=font_size))
  }

  if(!is.na(title)){
    stackedPlot <- stackedPlot +
      ggtitle(title)
    
    if(!is.na(font_size)){
      stackedPlot <- stackedPlot +
        theme(plot.title = element_text(size=font_size))
    }
  } 
  
  if(utils::packageVersion("ggplot2") >= '3.0.0'){
    stackedPlot <- stackedPlot +
      coord_flip(clip = "off")
  } else {
    stackedPlot <- stackedPlot +
      coord_flip()   
  }
  
  return(stackedPlot)
}