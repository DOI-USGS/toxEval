#' Grouped Boxplots
#' 
#' This function creates a set of boxplots based on the original input data 
#' modified by the processing steps above, and the choice of several input options. 
#' The result includes boxplots of the maximum or mean value per site based on 
#' the mean_logic argument This ensures that each site is represented equally 
#' regardless of how many samples are available per site. Boxplots are generated 
#' by the chosen category. Categories include "Biological", "Chemical Class", or 
#' "Chemical". Biological refers to the chosen toxCast annotation as defined 
#' in the groupCol argument of the \code{filter_groups} function. Chemical Class 
#' refers to the groupings of chemicals as defined in the "Class" column of the 
#' "Chemicals" tab of the input file. Choosing Chemical Class will generate a 
#' separate boxplot for each unique class. Choosing Chemical will generate a 
#' separate boxplot for each individual chemical in the data set.
#' 
#' The difference in the single-site output is that instead of listing the number 
#' of sites, it lists the number of unique chemical/endpoint combinations used to 
#' create the box plot.
#' 
#' @param chemicalSummary data frame from \code{get_chemical_summary}
#' @param category either "Biological", "Chemical Class", or "Chemical"
#' @param manual_remove vector of categories to remove
#' @param mean_logic logical.  TRUE takes the mean sample of each site,
#' FALSE takes the maximum sample of each site.
#' @param sum_logic logical. TRUE sums the EARs in a specified grouping,
#' FALSE does not. FALSE may be better for traditional benchmarks as
#' opposed to ToxCast benchmarks.
#' @param plot_ND logical whether or not to plot the non-detects
#' @param hit_threshold numeric threshold defining a "hit"
#' @param font_size numeric to adjust the axis font size
#' @param title character title for plot.
#' @param pallette vector of color pallette for fill. Can be a named vector
#' to specify specific color for specific categories. 
#' @export
#' @rdname plot_tox_boxplots
#' @import ggplot2
#' @importFrom stats median
#' @importFrom dplyr full_join filter mutate select left_join right_join
#' @examples
#' # This is the example workflow:
#' path_to_tox <-  system.file("extdata", package="toxEval")
#' file_name <- "OWC_data_fromSup.xlsx"
#'
#' full_path <- file.path(path_to_tox, file_name)
#' 
#' \dontrun{
#' tox_list <- create_toxEval(full_path)
#' ACClong <- get_ACC(tox_list$chem_info$CAS)
#' ACClong <- remove_flags(ACClong)
#' 
#' cleaned_ep <- clean_endPoint_info(endPointInfo)
#' filtered_ep <- filter_groups(cleaned_ep)
#' chemicalSummary <- get_chemical_summary(tox_list, ACClong, filtered_ep)
#' plot_tox_boxplots(chemicalSummary, "Biological")   
#' plot_tox_boxplots(chemicalSummary, "Chemical Class")
#' plot_tox_boxplots(chemicalSummary, "Chemical")
#' 
#'  # To turn off clipping:
#' class_plot <- plot_tox_boxplots(chemicalSummary, "Chemical Class")
#' gb <- ggplot2::ggplot_build(class_plot)
#' gt <- ggplot2::ggplot_gtable(gb)
#' 
#' gt$layout$clip[gt$layout$name=="panel"] <- "off"
#' 
#' grid::grid.draw(gt) 
#' 
#' cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", 
#'                "#0072B2", "#D55E00", "#CC79A7")
#' graphData <- tox_boxplot_data(chemicalSummary = chemicalSummary,
#'                               category = "Biological") 
#' cbValues <- colorRampPalette(cbPalette)(length(levels(graphData$category)))
#' names(cbValues) <- levels(graphData$category)
#' 
#' plot_tox_boxplots(chemicalSummary, 
#'                   hit_threshold = 0.1,
#'                   category = "Biological",
#'                   pallette = cbValues,
#'                   title = 'Maximum EAR per site, grouped by biological activity groupings') 
#' 
#' 
#' }
plot_tox_boxplots <- function(chemicalSummary, 
                              category = "Biological",
                              manual_remove = NULL,
                              mean_logic = FALSE,
                              sum_logic = TRUE,
                              plot_ND = TRUE, 
                              font_size = NA,
                              title = NA,
                              pallette = NA,
                              hit_threshold = NA){
  
  match.arg(category, c("Biological","Chemical Class","Chemical"))
  
  site <- EAR <- sumEAR <- meanEAR <- groupCol <- nonZero <- ".dplyr"
  x <- y <- CAS <- ".dplyr"

  if(category == "Chemical"){

    chemPlot <- plot_chemical_boxplots(chemicalSummary, 
                                       mean_logic = mean_logic,
                                       sum_logic = sum_logic,
                                       plot_ND = plot_ND,
                                       font_size = font_size,
                                       title = title,
                                       pallette = pallette,
                                       hit_threshold = hit_threshold)
    return(chemPlot)
    
  } else {
    
    if(!plot_ND){
      chemicalSummary <- chemicalSummary[chemicalSummary$EAR > 0,]
    }
    single_site <- length(unique(chemicalSummary$site)) == 1
    
    y_label <- fancyLabels(category, mean_logic, sum_logic, single_site)
    
    if(single_site){
      
      if(category == "Biological"){
        chemicalSummary$category <- chemicalSummary$Bio_category
      } else {
        chemicalSummary$category <- chemicalSummary$Class
      }
      
      pretty_logs_new <- prettyLogs(chemicalSummary$EAR)

      countNonZero <- chemicalSummary %>%
        group_by(category) %>%
        summarise(nonZero = as.character(length(unique(CAS))),
                  hits = as.character(sum(EAR > hit_threshold))) %>%
        data.frame() 
      
      countNonZero$hits[countNonZero$hits == "0"] <- ""
      
      label <- "# Chemicals"
      
      if(!is.null(manual_remove)){
        chemicalSummary <- filter(chemicalSummary, !(category %in% manual_remove))
      }
      
      orderColsBy <- chemicalSummary %>%
        group_by(category) %>%
        summarise(median = median(EAR[EAR != 0])) %>%
        arrange(median)
      
      orderedLevels <- orderColsBy$category
      
      if(any(is.na(orderColsBy$median))){
        orderedLevels <- c(orderColsBy$category[is.na(orderColsBy$median)],
                           orderColsBy$category[!is.na(orderColsBy$median)])
      }
      
      chemicalSummary$category <- factor(chemicalSummary$category,
                                         levels = orderedLevels[orderedLevels %in% chemicalSummary$category])
      
      bioPlot <- ggplot(data = chemicalSummary)+
        theme_bw() +
        xlab("") +
        theme(plot.background = element_rect(fill = "transparent",colour = NA),
              axis.text.y = element_text(color = "black", vjust = 0.2), 
              axis.text.x = element_text(color = "black", vjust = 0, margin = margin(-0.5,0,0,0)),
              panel.border = element_blank(),
              axis.ticks = element_blank(),
              plot.title = element_text(hjust = 0.5)) + 
        scale_y_log10(y_label,labels=fancyNumbers,breaks=pretty_logs_new) +
        geom_hline(yintercept = hit_threshold, linetype="dashed", color="black")
      
      if(!all(is.na(pallette))){
        bioPlot <- bioPlot +
          geom_boxplot(aes(x=category, y=EAR),lwd=0.1,outlier.size=1, fill = "steelblue") +
          scale_fill_manual(values = pallette) +
          theme(legend.position = "none")
      } else {
        bioPlot <- bioPlot +
          geom_boxplot(aes(x=category, y=EAR),lwd=0.1,outlier.size=1, fill = "steelblue") 
      }
      
    } else {
      
      graphData <- tox_boxplot_data(chemicalSummary = chemicalSummary,
                             category = category,
                             manual_remove = manual_remove,
                             mean_logic = mean_logic,
                             sum_logic = sum_logic)
      
      pretty_logs_new <- prettyLogs(graphData$meanEAR)
      
      countNonZero <- graphData %>%
        group_by(category) %>%
        summarise(nonZero = as.character(length(unique(site[meanEAR>0]))),
                  hits = as.character(sum(meanEAR > hit_threshold))) %>%
        data.frame()
      
      countNonZero$hits[countNonZero$hits == "0"] <- ""
      
      label <- "# Sites"
      
      bioPlot <- ggplot(data = graphData)+
        scale_y_log10(y_label,labels=fancyNumbers,breaks=pretty_logs_new) +
        theme_bw() +
        xlab("") +
        theme(plot.background = element_rect(fill = "transparent",colour = NA),
              axis.text.y = element_text(color = "black", vjust = 0.2), 
              axis.text.x = element_text(color = "black", vjust = 0),
              panel.border = element_blank(),
              axis.ticks = element_blank(),
              plot.title = element_text(hjust = 0.5, vjust = 0, margin = margin(-0.5,0,0,0))) +  
        geom_hline(yintercept = hit_threshold, linetype="dashed", color="black")
    
      if(!all(is.na(pallette))){
        bioPlot <- bioPlot +
          geom_boxplot(aes(x=category, y=meanEAR, fill = category),lwd=0.1,outlier.size=1) +
          scale_fill_manual(values = pallette) +
          theme(legend.position = "none")
      } else {
        bioPlot <- bioPlot +
          geom_boxplot(aes(x=category, y=meanEAR),lwd=0.1,outlier.size=1, fill = "steelblue")      
      }
    }
    if(!is.na(font_size)){
      bioPlot <- bioPlot +
        theme(axis.text = element_text(size = font_size),
              axis.title =   element_text(size=font_size))
    }
    
    if(packageVersion("ggplot2") >= '2.2.1.9000'){
      bioPlot <- bioPlot +
        coord_flip(clip = "off")
    } else {
      bioPlot <- bioPlot +
        coord_flip()      
    }
    
    plot_info <- ggplot_build(bioPlot)
    layout_stuff <- plot_info$layout
    
    if(packageVersion("ggplot2") >= "2.2.1.9000"){
      xmin <- 10^(layout_stuff$panel_scales_y[[1]]$range$range[1])
      xmax <- 10^(layout_stuff$panel_scales_y[[1]]$range$range[2])
      ymax <- layout_stuff$panel_scales_x[[1]]$range$range[2]
    } else {
      xmin <- suppressWarnings(10^(layout_stuff$panel_ranges[[1]]$x.range[1]))
      xmax <- suppressWarnings(10^(layout_stuff$panel_ranges[[1]]$x.range[2]))
      ymax <- suppressWarnings(layout_stuff$panel_ranges[[1]]$y.range[2])
    }
    
    bioPlot_w_labels <- bioPlot + 
      geom_text(data=countNonZero, aes(x=category, y=xmin,label=nonZero),size=ifelse(is.na(font_size),3,0.30*font_size)) +
      geom_text(data=data.frame(x = Inf, y=xmin, label = label, stringsAsFactors = FALSE), 
                aes(x = x,  y=y, label = label),
                size=ifelse(is.na(font_size),3,0.30*font_size))
      
    
    nHitsEP <- countNonZero$hits
    
    if(isTRUE(sum(as.numeric(nHitsEP), na.rm = TRUE) > 0)) {
      bioPlot_w_labels <- bioPlot_w_labels +
        geom_text(data=countNonZero, aes(x=category, y=ymax,label=nHitsEP),size=ifelse(is.na(font_size),3,0.30*font_size)) +
        geom_text(data=data.frame(x = Inf, y=ymax, label = "# Hits", stringsAsFactors = FALSE), 
                  aes(x = x,  y=y, label = label),
                  size=ifelse(is.na(font_size),3,0.30*font_size))
    }
    
    if(!is.na(hit_threshold)) {
      bioPlot_w_labels <- bioPlot_w_labels +
        geom_text(data=data.frame(x = Inf, y=hit_threshold, label = "Threshold", stringsAsFactors = FALSE), 
                  aes(x = x,  y=y, label = label),
                  size=ifelse(is.na(font_size),3,0.30*font_size))
    }
    
    if(!all(is.na(title))){
      bioPlot_w_labels <- bioPlot_w_labels +
        ggtitle(title) 
      
      if(!is.na(font_size)){
        bioPlot_w_labels <- bioPlot_w_labels +
          theme(plot.title = element_text(size=font_size, margin = margin(b = 5)))
      } else {
        bioPlot_w_labels <- bioPlot_w_labels +
          theme(plot.title = element_text(size=font_size, margin = margin(b = 5)))        
      }
    } 
    
    return(bioPlot_w_labels)
  }
  
}

#' @rdname plot_tox_boxplots
#' @export
tox_boxplot_data <- function(chemicalSummary, 
                      category = "Biological",
                      manual_remove = NULL, 
                      mean_logic = FALSE,
                      sum_logic = TRUE){
  
  match.arg(category, c("Biological","Chemical Class","Chemical"))

  if(category == "Chemical"){
    chm_sum <- graph_chem_data(chemicalSummary = chemicalSummary,
                               mean_logic = mean_logic,
                               sum_logic = sum_logic,
                               manual_remove = manual_remove)
    return(chm_sum)
  }
  
  site <- EAR <- sumEAR <- meanEAR <- groupCol <- nonZero <- ".dplyr"
  
  if(category == "Biological"){
    chemicalSummary$category <- chemicalSummary$Bio_category
  } else {
    chemicalSummary$category <- chemicalSummary$Class
  }

  if(!sum_logic){
    tox_boxplot_data <- chemicalSummary %>%
      group_by(site,category) %>%
      summarise(meanEAR=ifelse(mean_logic, mean(EAR), max(EAR))) %>%
      data.frame() 
      
  } else {
    tox_boxplot_data <- chemicalSummary %>%
      group_by(site,date,category) %>%
      summarise(sumEAR=sum(EAR)) %>%
      data.frame() %>%
      group_by(site, category) %>%
      summarise(meanEAR=ifelse(mean_logic, mean(sumEAR), max(sumEAR))) %>%
      data.frame()     
  }

  if(!is.null(manual_remove)){
    tox_boxplot_data <- filter(tox_boxplot_data, !(category %in% manual_remove))
  }
  
  orderColsBy <- tox_boxplot_data %>%
    group_by(category) %>%
    summarise(median = median(meanEAR[meanEAR != 0])) %>%
    arrange(median)
  
  orderedLevels <- orderColsBy$category
  
  if(any(is.na(orderColsBy$median))){
    orderedLevels <- c(orderColsBy$category[is.na(orderColsBy$median)],
                       orderColsBy$category[!is.na(orderColsBy$median)])
  }
  
  tox_boxplot_data$category <- factor(as.character(tox_boxplot_data$category), 
                               levels=orderedLevels[orderedLevels %in% unique(as.character(tox_boxplot_data$category))])
  
  return(tox_boxplot_data)
}

prettyLogs <- function(x){
  pretty_range <- range(x[x > 0])
  pretty_logs <- 10^(-10:10)
  log_index <- which(pretty_logs < pretty_range[2] & pretty_logs > pretty_range[1])
  log_index <- c(log_index[1]-1,log_index, log_index[length(log_index)]+1)
  pretty_logs_new <-  pretty_logs[log_index] 
  return(pretty_logs_new)
}


fancyNumbers2 <- function(n){
  textReturn <-  signif(n,digits = 2)
  textReturn <- as.character(textReturn)
  textReturn[length(textReturn)] <- paste(">",textReturn[length(textReturn)])
  textReturn[1] <- paste("<",textReturn[1])
  return(textReturn)
}

fancyNumbers <- function(n){
  nNoNA <- n[!is.na(n)]
  x <-gsub(pattern = "1e",replacement = "10^",x = format(nNoNA, scientific = TRUE))
  exponents <- as.numeric(sapply(strsplit(x, "\\^"), function(j) j[2]))

  base <- ifelse(exponents == 0, "1", ifelse(exponents == 1, "10","10^"))
  exponents[base == "1" | base == "10"] <- ""
  textNums <- rep(NA, length(n))  
  textNums[!is.na(n)] <- paste0(base,exponents)
  
  textReturn <- parse(text=textNums)
  return(textReturn)
}

fancyLabels <- function(category, mean_logic, sum_logic, single_site, sep=FALSE){
  
  pretty_cat <- switch(category, 
                       "Chemical" = "i = chemicals, j = samples, k = sites",
                       "Biological" = "i = chemicals in a specified grouping, j = samples, k = sites",
                       "Chemical Class" = "i = chemicals in a specified class, j = samples, k = sites"
  )
  
  word_stat <- ifelse(mean_logic, "mean", "max")

  if(single_site){
    
    y_label <- switch(category,
                      "Chemical Class" = "All EARs within a chemical class",
                      "Biological" = "All EARs within a grouping",
                      "Chemical" = "All EARs for each chemical")
    
    if(sep){
      y_label <- list(y_label=y_label, caption = NULL)
    }

    
  } else {
  
    if(sep){
      y_label <- bquote(.(word_stat) ~ 
                               group("[", 
                                     group("(",
                                           sum(" "  ~ EAR["[" *i* "]"]),
                                           ")")["[" *j* "]"],
                                     "]")
                             ["[" *k* "]"])
      y_label <- list(y_label = y_label, caption = pretty_cat)
    } else {

      y_label <- bquote(atop(.(word_stat) ~ 
                               group("[", 
                                     group("(",
                                           sum(" "  ~ EAR["[" *i* "]"]),
                                           ")")["[" *j* "]"],
                                     "]")
                             ["[" *k* "]"], .(pretty_cat)))
    }
  }
    
  return(y_label)
}
  