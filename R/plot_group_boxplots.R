#' Grouped Boxplots
#' 
#' The \code{plot_tox_boxplots} function creates a set of boxplots representing EAR 
#' values computed with the \code{\link{get_chemical_summary}} function, and 
#' dependent on the choice of several input options. See "Summarizing the data" 
#' in the Introduction vignette: \href{../doc/Introduction.html#summarize_data}{\code{vignette("Introduction", package = "toxEval")}}. 
#' for a description of how the EAR values are computed, aggregated,
#' and summarized. Choosing "Chemical Class" in the category argument 
#' will generate separate boxplots for each unique class. "Chemical" will generate 
#' boxplots for each individual chemical, and "Biological" will generate boxplots 
#' for each group in the selected ToxCast annotation.
#' 
#' It is also possible to display a threshold line using the hit_threshold argument. 
#' The graph will then include the number of sites with detections, the threshold 
#' line, and the number of "hits" indicating how many sites that have EAR values 
#' exceeding the hit_threshold.
#' 
#' The graph shows a slightly different result for a single site. For a single 
#' site graph, the number of chemicals that were detected and have associated endpoint 
#' ACCs represented are displayed.
#' 
#' The functions \code{plot_tox_boxplots} and \code{graph_chem_data} are functions that perform
#' the statistical calculations to create the plot. \code{graph_chem_data} is specific
#' to the "Chemical" plot, and \code{plot_tox_boxplots} is for "Biological" and 
#' "Chemical Class".
#' 
#' Box plots are standard Tukey representations. See "Box plot details" in the Basic Workflow vignette: 
#' \href{../doc/basicWorkflow.html#box_plot_details}{\code{vignette("basicWorkflow", package = "toxEval")}}
#' for more information.
#' 
#' @param chemical_summary Data frame from \code{\link{get_chemical_summary}}.
#' @param category Character. Either "Biological", "Chemical Class", or "Chemical".
#' @param manual_remove Vector of categories to remove.
#' @param mean_logic Logical.  \code{TRUE} displays the mean sample from each site,
#' \code{FALSE} displays the maximum sample from each site.
#' @param sum_logic Logical. \code{TRUE} sums the EARs in a specified grouping,
#' \code{FALSE} does not. \code{FALSE} may be better for traditional benchmarks as
#' opposed to ToxCast benchmarks.
#' @param plot_ND Logical. Logical. Whether or not to plot "Biological" groupings,
#' "Chemical Class" groupings, or "Chemical" that do not have any detections. 
#' @param hit_threshold Numeric threshold defining a "hit".
#' @param font_size Numeric value to adjust the axis font size.
#' @param title Character Title for plot.
#' @param palette Vector of color palette for boxplot fill. Can be a named vector
#' to specify specific colors for specific categories. 
#' @export
#' @rdname plot_tox_boxplots
#' @import ggplot2
#' @importFrom stats median
#' @examples
#' # This is the example workflow:
#' path_to_tox <-  system.file("extdata", package="toxEval")
#' file_name <- "OWC_data_fromSup.xlsx"
#'
#' full_path <- file.path(path_to_tox, file_name)
#' 
#' tox_list <- create_toxEval(full_path)
#' ACC <- get_ACC(tox_list$chem_info$CAS)
#' ACC <- remove_flags(ACC)
#' 
#' cleaned_ep <- clean_endPoint_info(end_point_info)
#' filtered_ep <- filter_groups(cleaned_ep)
#' chemical_summary <- get_chemical_summary(tox_list, ACC, filtered_ep)
#' plot_tox_boxplots(chemical_summary, "Biological")   
#' plot_tox_boxplots(chemical_summary, "Chemical Class")
#' plot_tox_boxplots(chemical_summary, "Chemical")
#' 
#' 
#' cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", 
#'                "#0072B2", "#D55E00", "#CC79A7")
#' graphData <- tox_boxplot_data(chemical_summary = chemical_summary,
#'                               category = "Biological") 
#' cbValues <- colorRampPalette(cbPalette)(length(levels(graphData$category)))
#' names(cbValues) <- levels(graphData$category)
#' 
#' plot_tox_boxplots(chemical_summary, 
#'                   hit_threshold = 0.1,
#'                   category = "Biological",
#'                   palette = cbValues,
#'                   title = 'Maximum EAR per site, grouped by biological activity groupings') 
#' 
plot_tox_boxplots <- function(chemical_summary, 
                              category = "Biological",
                              manual_remove = NULL,
                              mean_logic = FALSE,
                              sum_logic = TRUE,
                              plot_ND = TRUE, 
                              font_size = NA,
                              title = NA,
                              palette = NA,
                              hit_threshold = NA){
  
  match.arg(category, c("Biological","Chemical Class","Chemical"))
  
  site <- EAR <- sumEAR <- meanEAR <- groupCol <- nonZero <- ".dplyr"
  x <- y <- CAS <- ".dplyr"

  if(category == "Chemical"){

    chemPlot <- plot_chemical_boxplots(chemical_summary, 
                                       mean_logic = mean_logic,
                                       sum_logic = sum_logic,
                                       plot_ND = plot_ND,
                                       font_size = font_size,
                                       title = title,
                                       palette = palette,
                                       hit_threshold = hit_threshold)
    return(chemPlot)
    
  } else {
    
    if(!plot_ND){
      chemical_summary <- chemical_summary[chemical_summary$EAR > 0,]
    }
    single_site <- length(unique(chemical_summary$site)) == 1
    
    y_label <- fancyLabels(category, mean_logic, sum_logic, single_site)
    
    if(single_site){
      
      if(category == "Biological"){
        chemical_summary$category <- chemical_summary$Bio_category
      } else {
        chemical_summary$category <- chemical_summary$Class
      }
      
      pretty_logs_new <- prettyLogs(chemical_summary$EAR)

      countNonZero <- chemical_summary %>%
        dplyr::group_by(category) %>%
        dplyr::summarise(nonZero = as.character(length(unique(CAS))),
                  hits = as.character(sum(EAR > hit_threshold))) %>%
        data.frame() 
      
      countNonZero$hits[countNonZero$hits == "0"] <- ""
      
      label <- "# Chemicals"
      
      if(!is.null(manual_remove)){
        chemical_summary <- dplyr::filter(chemical_summary, !(category %in% manual_remove))
      }
      
      orderColsBy <- chemical_summary %>%
        dplyr::group_by(category) %>%
        dplyr::summarise(median = median(EAR[EAR != 0])) %>%
        dplyr::arrange(median)
      
      orderedLevels <- orderColsBy$category
      
      if(any(is.na(orderColsBy$median))){
        orderedLevels <- c(orderColsBy$category[is.na(orderColsBy$median)],
                           orderColsBy$category[!is.na(orderColsBy$median)])
      }
      
      chemical_summary$category <- factor(chemical_summary$category,
                                         levels = orderedLevels[orderedLevels %in% chemical_summary$category])
      
      bioPlot <- ggplot(data = chemical_summary)+
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
      
      if(!all(is.na(palette))){
        bioPlot <- bioPlot +
          geom_boxplot(aes(x=category, y=EAR),lwd=0.1,outlier.size=1, fill = "steelblue") +
          scale_fill_manual(values = palette) +
          theme(legend.position = "none")
      } else {
        bioPlot <- bioPlot +
          geom_boxplot(aes(x=category, y=EAR),lwd=0.1,outlier.size=1, fill = "steelblue") 
      }
      
    } else {
      
      graphData <- tox_boxplot_data(chemical_summary = chemical_summary,
                             category = category,
                             manual_remove = manual_remove,
                             mean_logic = mean_logic,
                             sum_logic = sum_logic)
      
      pretty_logs_new <- prettyLogs(graphData$meanEAR)
      
      countNonZero <- graphData %>%
        dplyr::group_by(category) %>%
        dplyr::summarise(nonZero = as.character(length(unique(site[meanEAR>0]))),
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
    
      if(!all(is.na(palette))){
        bioPlot <- bioPlot +
          geom_boxplot(aes(x=category, y=meanEAR, fill = category),lwd=0.1,outlier.size=1) +
          scale_fill_manual(values = palette) +
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
    
    if(utils::packageVersion("ggplot2") >= '3.0.0'){
      bioPlot <- bioPlot +
        coord_flip(clip = "off")
    } else {
      bioPlot <- bioPlot +
        coord_flip()      
    }
    
    plot_info <- ggplot_build(bioPlot)
    layout_stuff <- plot_info$layout
    
    if(utils::packageVersion("ggplot2") >= "3.0.0"){
      xmin <- 10^(layout_stuff$panel_scales_y[[1]]$range$range[1])
      xmax <- 10^(layout_stuff$panel_scales_y[[1]]$range$range[2])
      ymax <- length(layout_stuff$panel_scales_x[[1]]$range$range)
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
tox_boxplot_data <- function(chemical_summary, 
                      category = "Biological",
                      manual_remove = NULL, 
                      mean_logic = FALSE,
                      sum_logic = TRUE){
  
  match.arg(category, c("Biological","Chemical Class","Chemical"))

  if(category == "Chemical"){
    chm_sum <- graph_chem_data(chemical_summary = chemical_summary,
                               mean_logic = mean_logic,
                               sum_logic = sum_logic,
                               manual_remove = manual_remove)
    return(chm_sum)
  }
  
  site <- EAR <- sumEAR <- meanEAR <- groupCol <- nonZero <- ".dplyr"
  
  if(category == "Biological"){
    chemical_summary$category <- chemical_summary$Bio_category
  } else {
    chemical_summary$category <- chemical_summary$Class
  }

  if(!sum_logic){
    tox_boxplot_data <- chemical_summary %>%
      dplyr::group_by(site,category) %>%
      dplyr::summarise(meanEAR=ifelse(mean_logic, mean(EAR), max(EAR))) %>%
      data.frame() 
      
  } else {
    tox_boxplot_data <- chemical_summary %>%
      dplyr::group_by(site,date,category) %>%
      dplyr::summarise(sumEAR=sum(EAR)) %>%
      data.frame() %>%
      dplyr::group_by(site, category) %>%
      dplyr::summarise(meanEAR=ifelse(mean_logic, mean(sumEAR), max(sumEAR))) %>%
      data.frame()     
  }

  if(!is.null(manual_remove)){
    tox_boxplot_data <- dplyr::filter(tox_boxplot_data, !(category %in% manual_remove))
  }
  
  orderColsBy <- tox_boxplot_data %>%
    dplyr::group_by(category) %>%
    dplyr::summarise(median = median(meanEAR[meanEAR != 0])) %>%
    dplyr::arrange(median)
  
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
  # textReturn[length(textReturn)] <- paste(">",textReturn[length(textReturn)])
  # textReturn[1] <- paste("<",textReturn[1])
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

fancyLabels <- function(category, mean_logic, sum_logic, single_site, sep=FALSE, include_site = TRUE){
  
  pretty_cat <- switch(category, 
                       "Chemical" = "i = chemicals, j = samples, k = sites",
                       "Biological" = "i = chemicals in a specified grouping, j = samples, k = sites",
                       "Chemical Class" = "i = chemicals in a specified class, j = samples, k = sites"
  )
  if(!include_site){
    pretty_cat <- gsub(", k = sites","",pretty_cat)
  }
  pretty_cat <- bquote(italic(.(pretty_cat)))
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
      if(sum_logic){
        if(include_site){
          y_label <- bquote(italic(.(word_stat)) ~ 
                              group("[", 
                                    group("(",
                                          sum(" "*EAR["[" *i* "]"]),
                                          ")")["[" *j* "]"],
                                    "]")
                            ["[" *k* "]"])          
        } else {
          y_label <- bquote(italic(.(word_stat)) ~ 
                              group("[", sum(" "*EAR["[" *i* "]"]),"]")["[" *j* "]"])
        }

      } else {
        if(include_site){
          y_label <- bquote(italic(.(word_stat))* 
                              group("[", 
                                    italic(max) * group("(",EAR["[" *i* "]"],")")["[" *j* "]"],
                                    "]")
                            ["[" *k* "]"])           
        } else {
          y_label <- bquote(italic(.(word_stat)) * 
                              group("[", 
                                    italic(max) * group("(",EAR["[" *i* "]"],")")["[" *j* "]"],
                                    "]")) 
        }
       
      }
      y_label <- list(y_label = y_label, caption = pretty_cat)
      
    } else {
      
      if(sum_logic){
        if(include_site){
          y_label <- bquote(atop(italic(.(word_stat)) * 
                                   group("[", 
                                         group("(",
                                               sum(" "*EAR["[" *i* "]"]),
                                               ")")["[" *j* "]"],
                                         "]")
                                 ["[" *k* "]"], .(pretty_cat)))
        } else {
          y_label <- bquote(atop(italic(.(word_stat)) * 
                                   group("[", sum(" "*EAR["[" *i* "]"]),
                                               ")")["[" *j* "]"], 
                                 .(pretty_cat)))          
        }
      } else {

        if(include_site){
          y_label <- bquote(atop(italic(.(word_stat)) * 
                                   group("[", 
                                         italic(max) * group("(",EAR["[" *i* "]"],")")["[" *j* "]"],
                                         "]")
                                 ["[" *k* "]"], .(pretty_cat)))
        } else {
          y_label <- bquote(atop(italic(.(word_stat)) * 
                            italic(max) * group("(",EAR["[" *i* "]"],")")["[" *j* "]"], 
                            .(pretty_cat)))
        }
        
      }
    }
  }
    
  return(y_label)
}
  