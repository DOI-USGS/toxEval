#' plot_tox_boxplots
#' 
#' Plot boxplot of groups
#' @param chemicalSummary data frame from \code{get_chemical_summary}
#' @param category either "Biological", "Chemical Class", or "Chemical"
#' @param manual_remove vector of categories to remove
#' @param mean_logic logical \code{TRUE} is mean, \code{FALSE} is maximum
#' @param plot_ND logical whether or not to plot the non-detects
#' @param font_size numeric to adjust the axis font size
#' @param title character title for plot. 
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
#' }
plot_tox_boxplots <- function(chemicalSummary, 
                              category = "Biological",
                              manual_remove = NULL,
                              mean_logic = FALSE,
                              plot_ND = TRUE, 
                              font_size = NA,
                              title = NA){
  
  match.arg(category, c("Biological","Chemical Class","Chemical"))

  site <- EAR <- sumEAR <- meanEAR <- groupCol <- nonZero <- ".dplyr"
  x <- y <- CAS <- ".dplyr"

  if(category == "Chemical"){

    chemPlot <- plot_chemical_boxplots(chemicalSummary, 
                                       mean_logic = mean_logic, 
                                       plot_ND = plot_ND,
                                       font_size = font_size,
                                       title = title)
    return(chemPlot)
    
  } else {
    
    if(!plot_ND){
      chemicalSummary <- chemicalSummary[chemicalSummary$EAR > 0,]
    }
    single_site <- length(unique(chemicalSummary$site)) == 1
    
    if(single_site){
      
      if(category == "Biological"){
        chemicalSummary$category <- chemicalSummary$Bio_category
      } else {
        chemicalSummary$category <- chemicalSummary$Class
      }
      
      countNonZero <- chemicalSummary %>%
        group_by(category) %>%
        summarise(nonZero = as.character(length(unique(CAS)))) %>%
        data.frame() 
      
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
        coord_flip() +
        theme_bw() +
        xlab("") +
        theme(plot.background = element_rect(fill = "transparent",colour = NA),
              axis.text.y = element_text(color = "black", vjust = 0.2), 
              axis.text.x = element_text(color = "black", vjust = 0, margin = margin(-0.5,0,0,0)),
              panel.border = element_blank(),
              axis.ticks = element_blank()) + 
        geom_boxplot(aes(x=category, y=EAR),lwd=0.1,outlier.size=1, fill = "steelblue") +
        scale_y_log10("EAR Per Sample",labels=fancyNumbers) 
      
    } else {
      graphData <- tox_boxplot_data(chemicalSummary = chemicalSummary,
                             category = category,
                             manual_remove = manual_remove,
                             mean_logic = mean_logic)
      
      countNonZero <- graphData %>%
        group_by(category) %>%
        summarise(nonZero = as.character(length(unique(site[meanEAR>0])))) %>%
        data.frame() 
      
      label <- "# Sites"
      
      bioPlot <- ggplot(data = graphData)+
        coord_flip() +
        theme_bw() +
        xlab("") +
        theme(plot.background = element_rect(fill = "transparent",colour = NA),
              axis.text.y = element_text(color = "black", vjust = 0.2), 
              axis.text.x = element_text(color = "black", vjust = 0, margin = margin(-0.5,0,0,0)),
              panel.border = element_blank(),
              axis.ticks = element_blank()) +  
        geom_boxplot(aes(x=category, y=meanEAR),lwd=0.1,outlier.size=1, fill = "steelblue") +
        scale_y_log10("Maximum EAR Per Site",labels=fancyNumbers) 
    }
    
    if(!is.na(font_size)){
      bioPlot <- bioPlot +
        theme(axis.text = element_text(size = font_size))
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
    
    if(!is.na(title)){
      bioPlot_w_labels <- bioPlot_w_labels +
        ggtitle(title)
      
      if(!is.na(font_size)){
        bioPlot_w_labels <- bioPlot_w_labels +
          theme(plot.title = element_text(size=font_size))
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
                      mean_logic = FALSE){
  
  match.arg(category, c("Biological","Chemical Class"))
  
  site <- EAR <- sumEAR <- meanEAR <- groupCol <- nonZero <- ".dplyr"

  
  if(category == "Biological"){
    chemicalSummary$category <- chemicalSummary$Bio_category
  } else {
    chemicalSummary$category <- chemicalSummary$Class
  }
  
  tox_boxplot_data <- chemicalSummary %>%
    group_by(site,date,category) %>%
    summarise(sumEAR=sum(EAR)) %>%
    data.frame() %>%
    group_by(site, category) %>%
    summarise(meanEAR=ifelse(mean_logic, mean(sumEAR), max(sumEAR))) %>%
    data.frame() 
  
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

  