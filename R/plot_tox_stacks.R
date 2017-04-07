#' plot_tox_boxplots
#' 
#' Plot boxplot of groups
#' @param chemicalSummary data frame from \code{get_chemical_summary}
#' @param category either "Biological", "Chemical Class", or "Chemical"
#' @param chem_site data frame with at least columns SiteID, site_grouping,  and Short Name
#' @param mean_logic logical \code{TRUE} is mean, \code{FALSE} is maximum
#' @param manual_remove vector of categories to remove
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
#' plot_tox_stacks(chemicalSummary, chem_site, "Biological")   
#' plot_tox_stacks(chemicalSummary, chem_site, "Chemical Class")
#' plot_tox_stacks(chemicalSummary, chem_site, "Chemical") 
plot_tox_stacks <- function(chemicalSummary, 
                            chem_site,
                            category = "Biological",
                            mean_logic = FALSE,
                            manual_remove = NULL){
  
  match.arg(category, c("Biological","Chemical Class","Chemical"))
  
  site <- EAR <- sumEAR <- meanEAR <- groupCol <- nonZero <- ".dplyr"
  
  graphData <- graphData(chemicalSummary = chemicalSummary,
                         category = category,
                         manual_remove = manual_remove)

  siteToFind <- unique(chemicalSummary$shortName)

  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  cbValues <- colorRampPalette(cbPalette)(length(levels(graphData$category)))
  set.seed(4)
  cbValues <- sample(cbValues)

  siteLimits <- chem_site$`Short Name`

  if(length(siteToFind) > 1){
    graphData <- graphData %>%
      left_join(select(chem_site, SiteID, site_grouping, `Short Name`),
                by=c("site"="SiteID"))

    upperPlot <- ggplot(graphData, 
                        aes(x=`Short Name`, y=meanEAR, fill = category)) +
      geom_bar(stat="identity") +
      theme_minimal() +
      theme(legend.title = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1)) +
      xlab("") +
      ylab(paste(ifelse(mean_logic,"Mean","Maximum"), "EAR Per Site")) +
      scale_fill_manual(values = cbValues, drop=FALSE) +
      facet_grid(. ~ site_grouping, scales="free", space="free") 

  } else {

    chemGroupBPOneSite <- chemicalSummary %>%
      select(-site)

    upperPlot <- ggplot(chemGroupBPOneSite, aes(x=date, y=EAR, fill = category)) +
      geom_bar(stat="identity")+
      theme_minimal() +
      theme(axis.text.x=element_blank(),
            axis.ticks=element_blank(),
            legend.title = element_blank())+
      xlab("Individual Samples") +
      ylab("EAR") +
      scale_fill_discrete("", drop=FALSE) +
      scale_fill_manual(values = cbValues, drop=FALSE)

  }
  
  return(upperPlot)
}