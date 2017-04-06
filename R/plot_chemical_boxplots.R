#' plot_chemical_boxplots
#' 
#' Plot boxplot of chemicals
#' @param graphData data frame from \code{graph_chem_data}
#' @export
#' @import ggplot2
#' @importFrom stats median quantile
#' @importFrom dplyr full_join filter mutate select left_join right_join
#' @examples
#' library(readxl)
#' path_to_tox <-  system.file("extdata", package="toxEval")
#' file_name <- "OWC_data_fromSup.xlsx"
#' full_path <- file.path(path_to_tox, file_name)
#' 
#' chem_data <- read_excel(full_path, sheet = "Data")
#' chem_info <- read_excel(full_path, sheet = "Chemicals") 
#' chem_site <- read_excel(full_path, sheet = "Sites")
#' ACClong <- get_ACC(chem_info$CAS)
#' ACClong <- remove_flags(ACClong)
#' 
#' cleaned_ep <- clean_endPoint_info(endPointInfo)
#' filtered_ep <- filter_groups(cleaned_ep)
#' 
#' chemicalSummary <- get_chemical_summary(ACClong,
#'                                         filtered_ep,
#'                                        chem_data, 
#'                                         chem_site, 
#'                                         chem_info)
#' 
#' graphData <- graph_chem_data(chemicalSummary)
#' plot_chemical_boxplots(graphData)
plot_chemical_boxplots <- function(graphData){
  
  site <- EAR <- sumEAR <- meanEAR <- groupCol <- nonZero <- ".dplyr"
  chnm <- Class <- maxEAR <- ".dplyr"
    
  cbValues <- c("#DCDA4B","#999999","#00FFFF","#CEA226","#CC79A7","#4E26CE",
                "#FFFF00","#78C15A","#79AEAE","#FF0000","#00FF00","#B1611D",
                "#FFA500","#F4426e")

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

#' graph_chem_data
#' 
#' Get chemical data summarized for plots
#' @param chemicalSummary data frame from \code{graph_chem_data}
#' @param manual_remove vector of categories to remove
#' @export
#' @importFrom stats median
#' @importFrom dplyr full_join filter mutate select left_join right_join
#' @examples
#' library(readxl)
#' path_to_tox <-  system.file("extdata", package="toxEval")
#' file_name <- "OWC_data_fromSup.xlsx"
#' full_path <- file.path(path_to_tox, file_name)
#' 
#' chem_data <- read_excel(full_path, sheet = "Data")
#' chem_info <- read_excel(full_path, sheet = "Chemicals") 
#' chem_site <- read_excel(full_path, sheet = "Sites")
#' ACClong <- get_ACC(chem_info$CAS)
#' ACClong <- remove_flags(ACClong)
#' 
#' cleaned_ep <- clean_endPoint_info(endPointInfo)
#' filtered_ep <- filter_groups(cleaned_ep)
#' 
#' chemicalSummary <- get_chemical_summary(ACClong,
#'                                         filtered_ep,
#'                                        chem_data, 
#'                                         chem_site, 
#'                                         chem_info)
#' 
#' graphData <- graph_chem_data(chemicalSummary)
graph_chem_data <- function(chemicalSummary, manual_remove=NULL){
  
  site <- chnm <- Class <- EAR <- sumEAR <- maxEAR <- ".dplyr"
  
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
  
  return(graphData)
}

