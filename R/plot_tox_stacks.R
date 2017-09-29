#' plot_tox_boxplots
#' 
#' Plot boxplot of groups
#' @param chemicalSummary data frame from \code{get_chemical_summary}
#' @param category either "Biological", "Chemical Class", or "Chemical"
#' @param chem_site data frame with at least columns SiteID, site_grouping,  and Short Name
#' @param mean_logic logical \code{TRUE} is mean, \code{FALSE} is maximum
#' @param manual_remove vector of categories to remove
#' @param include_legend logical to include legend or not
#' @export
#' @import ggplot2
#' @importFrom stats median
#' @importFrom grDevices colorRampPalette
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
#' plot_tox_stacks(chemicalSummary, chem_site, "Biological")   
#' plot_tox_stacks(chemicalSummary, chem_site, "Chemical Class")
#' plot_tox_stacks(chemicalSummary, chem_site, "Chemical", include_legend = FALSE) 
plot_tox_stacks <- function(chemicalSummary, 
                            chem_site,
                            category = "Biological",
                            mean_logic = FALSE,
                            manual_remove = NULL,
                            include_legend = TRUE){
  
  match.arg(category, c("Biological","Chemical Class","Chemical"))
  
  site <- EAR <- sumEAR <- meanEAR <- groupCol <- nonZero <- ".dplyr"
  SiteID <- site_grouping <- `Short Name` <- ".dplyr"
    
  if(include_legend & category == "Chemical"){
    graphData <- graph_chem_data(chemicalSummary = chemicalSummary,
                           manual_remove = manual_remove,
                           mean_logic = mean_logic)   
    names(graphData)[names(graphData) == "maxEAR"] <- "meanEAR"
    names(graphData)[names(graphData) == "chnm"] <- "category"
  } else {
    graphData <- graphData(chemicalSummary = chemicalSummary,
                           category = category,
                           manual_remove = manual_remove,
                           mean_logic = mean_logic) 
    if(category == "Chemical"){
      graphData$category <- graphData$chnm
    } 
  }

  siteToFind <- unique(chemicalSummary$shortName)

  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  cbValues <- colorRampPalette(cbPalette)(length(levels(graphData$category)))
  set.seed(4)
  cbValues <- sample(cbValues)

  siteLimits <- chem_site$`Short Name`

  if(length(siteToFind) > 1){
    graphData <- graphData %>%
      left_join(chem_site[, c("SiteID", "site_grouping", "Short Name")],
                by=c("site"="SiteID"))

    upperPlot <- ggplot(graphData, 
                        aes(x=`Short Name`, y=meanEAR, fill = category)) +
      theme_minimal() +
      xlab("") +
      ylab(paste(ifelse(mean_logic,"Mean","Maximum"), "EAR Per Site")) +
      facet_grid(. ~ site_grouping, scales="free", space="free") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

  } else {

    graphData <- chemicalSummary %>%
      select(-site)

    if(category == "Biological"){
      graphData$category <- graphData$Bio_category
    } else if (category == "Chemical Class"){
      graphData$category <- graphData$Class
    } else {
      graphData$category <- graphData$chnm
    }
    
    upperPlot <- ggplot(graphData, aes(x=date, y=EAR, fill = category)) +
      theme_minimal() +
      theme(axis.text.x=element_blank(),
            axis.ticks.x=element_blank()) +
      xlab("Individual Samples") +
      ylab("EAR") 

  }
  
  upperPlot <- upperPlot +
    geom_bar(stat="identity") 
  
  if(include_legend && length(unique(graphData$category)) <= length(cbValues)){
    upperPlot <- upperPlot + 
      scale_fill_manual(name = category,values = cbValues, drop=FALSE)

  } else {
    upperPlot <- upperPlot +
      theme(legend.position="none")
  }
  
  return(upperPlot)
}