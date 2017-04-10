#' plot_tox_boxplots
#' 
#' Plot boxplot of groups
#' @param chemicalSummary data frame from \code{get_chemical_summary}
#' @param category either "Biological", "Chemical Class", or "Chemical"
#' @param filterBy character either "All" or one of the filtered categories.
#' @param manual_remove vector of categories to remove
#' @param mean_logic logical \code{TRUE} is mean, \code{FALSE} is maximum
#' @param hit_threshold numeric threshold defining a "hit"
#' @export
#' @import ggplot2
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
#'  plot_tox_endpoints(chemicalSummary, filterBy = "Cell Cycle")
plot_tox_endpoints <- function(chemicalSummary, 
                              category = "Biological",
                              filterBy = "All",
                              manual_remove = NULL,
                              hit_threshold = 0.1,
                              mean_logic = FALSE){
  
  match.arg(category, c("Biological","Chemical Class","Chemical"))
  
  site <- endPoint <- EAR <- sumEAR <- meanEAR <-  ".dplyr"
  
  if(category == "Biological"){
    chemicalSummary$category <- chemicalSummary$Bio_category
  } else if(category == "Chemical Class") {
    chemicalSummary$category <- chemicalSummary$Class
  } else {
    chemicalSummary$category <- chemicalSummary$chnm
  }
      
  graphData <- chemicalSummary %>%
      group_by(site,date,category,endPoint) %>%
      summarise(sumEAR=sum(EAR)) %>%
      data.frame() %>%
      group_by(site, category,endPoint) %>%
      summarise(meanEAR=ifelse(mean_logic,mean(sumEAR),max(sumEAR))) %>%
      data.frame() %>%
      mutate(category=as.character(category))
  
  if(filterBy != "All"){
    graphData <- graphData %>%
      filter_(paste0("category == '", filterBy,"'"))

    countNonZero <- graphData %>%
      group_by(endPoint) %>%
      summarise(nonZero = as.character(sum(meanEAR>0)),
                hits = as.character(sum(meanEAR > hit_threshold)))

    countNonZero$hits[countNonZero$hits == "0"] <- ""

    namesToPlotEP <- as.character(countNonZero$endPoint)
    nSamplesEP <- countNonZero$nonZero
    nHitsEP <- countNonZero$hits
  }

  orderColsBy <- graphData %>%
    group_by(endPoint) %>%
    summarise(median = quantile(meanEAR[meanEAR != 0],0.5)) %>%
    arrange(median)

  orderedLevelsEP <- orderColsBy$endPoint

  if(any(is.na(orderColsBy$median))){
    orderedLevelsEP <- c(orderColsBy$endPoint[is.na(orderColsBy$median)],
                        orderColsBy$endPoint[!is.na(orderColsBy$median)])
  }

  graphData$endPoint <- factor(graphData$endPoint, levels = orderedLevelsEP)

  stackedPlot <- ggplot(graphData)+
    scale_y_log10(paste(ifelse(mean_logic,"Mean","Maximum"), "EAR Per Site"),labels=fancyNumbers) +
    geom_boxplot(aes(x=endPoint, y=meanEAR), fill = "steelblue") +
    theme_minimal() +
    xlab("") +
    theme(axis.text.y = element_text(vjust = .25,hjust=1)) +
    geom_hline(yintercept = hit_threshold, linetype="dashed", color="black")

  if(filterBy != "All"){

    ymin <- 10^(ggplot_build(stackedPlot)$layout$panel_ranges[[1]]$y.range)[1]
    ymax <- 10^(ggplot_build(stackedPlot)$layout$panel_ranges[[1]]$y.range)[2]

    xmax <- ggplot_build(stackedPlot)$layout$panel_ranges[[1]]$x.range[2]
    xmin <- ggplot_build(stackedPlot)$layout$panel_ranges[[1]]$x.range[1]

    stackedPlot <- stackedPlot +
      geom_text(data=data.frame(), aes(x=namesToPlotEP, y=ymin,label=nSamplesEP),size=5) +
      geom_text(data=data.frame(), aes(x=namesToPlotEP, y=ymax,label=nHitsEP),size=5)

    df1 <- data.frame(y = c(ymin,hit_threshold,ymax), 
                      text = c("# Non Zero",
                               "Hit Threshold",
                               "# Hits"), stringsAsFactors = FALSE)

    for(i in 1:3){
      stackedPlot <- stackedPlot +
        annotation_custom(
          grob = textGrob(label = df1$text[i], gp = gpar(cex = 0.75)),
          ymin = log10(df1$y[i]),      # Vertical position of the textGrob
          ymax = log10(df1$y[i]),
          xmin = xmax+0.05,  # Note: The grobs are positioned outside the plot area
          xmax = xmax+0.05)
    }
  }
      
  stackedPlot <- stackedPlot +
        coord_flip()
  
  return(stackedPlot)
}