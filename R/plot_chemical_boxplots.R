#' plot_chemical_boxplots
#' 
#' Plot boxplot of chemicals
#' @param chemicalSummary data frame from \code{get_chemical_summary}
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
plot_chemical_boxplots <- function(chemicalSummary, 
                                filtered_ep,
                                manual_remove = NULL){
  
  site <- EAR <- sumEAR <- meanEAR <- groupCol <- nonZero <- ".dplyr"
  
  cbValues <- c("#DCDA4B","#999999","#00FFFF","#CEA226","#CC79A7","#4E26CE",
                "#FFFF00","#78C15A","#79AEAE","#FF0000","#00FF00","#B1611D",
                "#FFA500","#F4426e")
  
  graphData <- chemicalSummary %>%
    group_by(site,date,chnm, Class) %>%
    summarise(sumEAR=sum(EAR)) %>%
    data.frame() %>%
    group_by(site, chnm, Class) %>%
    summarise(maxEAR=max(sumEAR)) %>%
    data.frame() 
  
  if(!is.null(manual_remove)){
    graphData <- filter(graphData, !(chnm %in% manual_remove))
  }
  
  orderClass <- graphData %>%
    group_by(Class,chnm) %>%
    summarise(median = median(maxEAR[maxEAR != 0])) %>%
    data.frame() %>%
    arrange(desc(median)) %>%
    filter(!duplicated(Class)) %>%
    arrange(median)
  
  orderChem <- graphData %>%
    group_by(chnm,Class) %>%
    summarise(median = quantile(maxEAR[maxEAR != 0],0.5)) %>%
    data.frame() %>%
    mutate(Class = factor(Class, levels=orderClass$Class)) %>%
    arrange(Class, !is.na(median), median)
  
  orderedLevels <- as.character(orderChem$chnm)
  
  graphData <- graphData %>%
    mutate(chnm = factor(chnm, levels=orderedLevels)) %>%
    mutate(Class = factor(Class, levels = rev(levels(orderChem$Class))))
  
  countNonZero <- graphData %>%
    select(chnm, Class, maxEAR) %>%
    group_by(chnm, Class) %>%
    summarize(nonZero = as.character(sum(maxEAR>0)))
  
  toxPlot_All <- ggplot(data=graphData) +
    scale_y_log10(labels=fancyNumbers)  +
    geom_boxplot(aes(x=chnm, y=maxEAR, fill=Class),
                 lwd=0.1,outlier.size=1)  +
    theme_bw() +
    scale_x_discrete(drop=TRUE) +
    coord_flip() +
    theme(axis.text = element_text( color = "black"),
          axis.text.y = element_text(size=7),
          axis.title=element_blank(),
          panel.background = element_blank(),
          plot.background = element_rect(fill = "transparent",colour = NA),
          strip.background = element_rect(fill = "transparent",colour = NA),
          strip.text.y = element_blank()) +
    guides(fill=guide_legend(ncol=6)) +
    theme(legend.position="bottom",
          legend.justification = "left",
          legend.background = element_rect(fill = "transparent", colour = "transparent"),
          legend.title=element_blank(),
          legend.text = element_text(size=8),
          legend.key.height = unit(1,"line")) +
    scale_fill_manual(values = cbValues, drop=FALSE) 
  
  ymin <- 0.3*10^-6
  ymax <- ggplot_build(toxPlot_All)$layout$panel_ranges[[1]]$y.range[2]
  
  toxPlot_All_withLabels <- toxPlot_All +
    geom_text(data=countNonZero, aes(x=chnm,label=nonZero, y=ymin), size=2.5) 
  
  return(toxPlot_All_withLabels)
  
}




