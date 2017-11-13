#' plot_chemical_boxplots
#' 
#' Plot boxplot of chemicals
#' @param chemicalSummary data frame from \code{graph_chem_data}
#' @param manual_remove vector of categories to remove
#' @param mean_logic logical \code{TRUE} is mean, \code{FALSE} is maximum
#' @param plot_ND logical whether or not to plot the non-detects
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
#' exclusion <- read_excel(full_path, sheet = "Exclude")
#' 
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
#'                                         chem_info,
#'                                         exclusion)
#'                                         
#' plot_chemical_boxplots(chemicalSummary)
plot_chemical_boxplots <- function(chemicalSummary, 
                                   manual_remove=NULL,
                                   mean_logic = FALSE,
                                   plot_ND = TRUE){
  
  site <- EAR <- sumEAR <- meanEAR <- groupCol <- nonZero <- ".dplyr"
  chnm <- Class <- maxEAR <- ".dplyr"
    
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
  
  if(length(unique(chemicalSummary$site)) == 1){
    
    countNonZero <- chemicalSummary %>%
      select(chnm, Class, EAR) %>%
      group_by(chnm, Class) %>%
      summarize(nonZero = as.character(sum(EAR>0)))
    
    toxPlot_All <- ggplot(data=chemicalSummary) +
      geom_boxplot(aes(x=chnm, y=EAR, fill=Class),
                   lwd=0.1,outlier.size=1)  
    
  } else {
    graphData <- graph_chem_data(chemicalSummary=chemicalSummary, 
                                 manual_remove=manual_remove,
                                 mean_logic=mean_logic)
    
    countNonZero <- graphData %>%
      select(chnm, Class, maxEAR) %>%
      group_by(chnm, Class) %>%
      summarize(nonZero = as.character(sum(maxEAR>0)))
    
    toxPlot_All <- ggplot(data=graphData) +
      geom_boxplot(aes(x=chnm, y=maxEAR, fill=Class),
                   lwd=0.1,outlier.size=1)  
  }
  
  toxPlot_All <- toxPlot_All +
    scale_y_log10(labels=fancyNumbers)  +
    theme_bw() +
    scale_x_discrete(drop = TRUE) +
    coord_flip() +
    theme(axis.text = element_text( color = "black"),
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
          legend.key.height = unit(1,"line")) +
    scale_fill_manual(values = cbValues, drop=FALSE)

  plot_info <- ggplot_build(toxPlot_All)
  layout_stuff <- plot_info$layout
  
  if(packageVersion("ggplot2") >= "2.2.1.9000"){
    ymin <- 10^(layout_stuff$panel_scales_y[[1]]$range$range[1])
  } else {
    ymin <- 10^(layout_stuff$panel_ranges[[1]]$x.range[1])
  }
  
  toxPlot_All_withLabels <- toxPlot_All +
    geom_text(data=countNonZero, aes(x=chnm,label=nonZero, y=ymin), size=2.5)

  return(toxPlot_All_withLabels)
  
}

#' graph_chem_data
#' 
#' Get chemical data summarized for plots
#' @param chemicalSummary data frame from \code{graph_chem_data}
#' @param manual_remove vector of categories to remove
#' @param mean_logic logical \code{TRUE} is mean, \code{FALSE} is maximum
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
graph_chem_data <- function(chemicalSummary, 
                            manual_remove=NULL,
                            mean_logic = FALSE){
  
  site <- chnm <- Class <- EAR <- sumEAR <- maxEAR <- ".dplyr"

  graphData <- chemicalSummary %>%
    group_by(site,date,chnm, Class) %>%
    summarise(sumEAR=sum(EAR)) %>%
    data.frame() %>%
    group_by(site, chnm, Class) %>%
    summarise(maxEAR=ifelse(mean_logic,mean(sumEAR),max(sumEAR))) %>%
    data.frame() 
  
  if(!is.null(manual_remove)){
    graphData <- filter(graphData, !(chnm %in% manual_remove))
  }
  
  return(graphData)
}

