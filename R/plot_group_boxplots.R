#' plot_group_boxplots
#' 
#' Plot boxplot of groups
#' @param chemicalSummary data frame from \code{get_chemical_summary}
#' @param category either "Biological" or "Chemical Class"
#' @param filtered_ep data frame from \code{filter_groups}
#' @param manual_remove vector of categories to remove
#' @export
#' @import ggplot2
#' @importFrom stats median
#' @importFrom dplyr full_join filter mutate select left_join right_join
#' @examples
#' CAS <- c("121-00-6","136-85-6","80-05-7","84-65-1","5436-43-1","126-73-8")
#' ACClong <- get_ACC(CAS)
#' ACClong <- remove_flags(ACClong)
#' cleaned_ep <- clean_endPoint_info(endPointInfo)
#' filtered_ep <- filter_groups(cleaned_ep)
#' 
plot_group_boxplots <- function(chemicalSummary, 
                                filtered_ep,
                                category = "Biological",
                                manual_remove = NULL){
  
  match.arg(category, c("Biological","Chemical Class"))
  
  site <- EAR <- sumEAR <- meanEAR <- groupCol <- nonZero <- ".dplyr"
  
  if(category == "Biological"){
    chemicalSummary$category <- chemicalSummary$Bio_category
  } else {
    chemicalSummary$category <- chemicalSummary$Class
  }
  
  graphData <- chemicalSummary %>%
    group_by(site,date,category) %>%
    summarise(sumEAR=sum(EAR)) %>%
    data.frame() %>%
    group_by(site, category) %>%
    summarise(meanEAR=max(sumEAR)) %>%
    data.frame() 
  
  if(!is.null(manual_remove)){
    graphData <- filter(graphData, !(category %in% manual_remove))
  }
    
  
  orderColsBy <- graphData %>%
    group_by(category) %>%
    summarise(median = median(meanEAR[meanEAR != 0])) %>%
    arrange(median)
  
  orderedLevels <- orderColsBy$category
  
  if(any(is.na(orderColsBy$median))){
    orderedLevels <- c(orderColsBy$category[is.na(orderColsBy$median)],
                       orderColsBy$category[!is.na(orderColsBy$median)])
  }
  
  orderNames <- names(table(select(filtered_ep, groupCol)))
  
  orderedLevels <- c(orderNames[!(orderNames %in% orderedLevels)],orderedLevels)
  
  countNonZero <- graphData %>%
    group_by(category) %>%
    summarise(nonZero = as.character(length(unique(site[meanEAR>0])))) %>%
    data.frame() 
  
  
  
  graphData$category <- factor(as.character(graphData$category), levels=orderedLevels)
  
  bioPlot <- ggplot(graphData)+
    scale_y_log10("Maximum EAR Per Site",labels=fancyNumbers)+
    geom_boxplot(aes(x=category, y=meanEAR),lwd=0.1,outlier.size=1, fill = "steelblue") +
    coord_flip() +
    theme_bw() +
    xlab("") +
    theme(plot.background = element_rect(fill = "transparent",colour = NA),
          axis.text.y = element_text(size=10, color = "black", vjust = 0.2), 
          axis.text.x = element_text(size=10, color = "black", vjust = 0, margin = margin(-0.5,0,0,0)),
          axis.title = element_text(size=10))
  
  xmin <- suppressWarnings(10^(ggplot_build(bioPlot)$layout$panel_ranges[[1]]$x.range[1]))
  xmax <- suppressWarnings(10^(ggplot_build(bioPlot)$layout$panel_ranges[[1]]$x.range[2]))
  ymax <- suppressWarnings(ggplot_build(bioPlot)$layout$panel_ranges[[1]]$y.range[1])
  
  bioPlot <- bioPlot + 
    geom_text(data=countNonZero, aes(x=category, y=xmin,label=nonZero),size=3) 
  
  return(bioPlot)
}


fancyNumbers2 <- function(n){
  textReturn <-  signif(n,digits = 2)
  textReturn <- as.character(textReturn)
  textReturn[length(textReturn)] <- paste(">",textReturn[length(textReturn)])
  textReturn[1] <- paste("<",textReturn[1])
  return(textReturn)
}

#' fancyNumbers
#' 
#' Plot fancyNumbers of groups
#' @param n vector
#' @export
fancyNumbers <- function(n){
  nNoNA <- n[!is.na(n)]
  x <-gsub(pattern = "1e",replacement = "10^",x = format(nNoNA, scientific = TRUE))
  exponents <- as.numeric(sapply(strsplit(x, "\\^"), function(j) j[2]))
  # browser()
  base <- ifelse(exponents == 0, "1", ifelse(exponents == 1, "10","10^"))
  exponents[base == "1" | base == "10"] <- ""
  textNums <- rep(NA, length(n))  
  textNums[!is.na(n)] <- paste0(base,exponents)
  
  textReturn <- parse(text=textNums)
  return(textReturn)
}

  