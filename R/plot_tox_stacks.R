#' plot_tox_boxplots
#' 
#' Plot boxplot of groups
#' @param chemicalSummary data frame from \code{get_chemical_summary}
#' @param category either "Biological", "Chemical Class", or "Chemical"
#' @param chem_site data frame with at least columns SiteID, site_grouping,  and Short Name
#' @param mean_logic logical \code{TRUE} is mean, \code{FALSE} is maximum
#' @param manual_remove vector of categories to remove
#' @param include_legend logical to include legend or not
#' @param font_size numeric to adjust the axis font size
#' @param title character title for plot. 
#' @export
#' @import ggplot2
#' @importFrom stats median
#' @importFrom grDevices colorRampPalette
#' @importFrom dplyr full_join filter mutate select left_join right_join
#' @examples
#' # This is the example workflow:
#' path_to_tox <-  system.file("extdata", package="toxEval")
#' file_name <- "OWC_data_fromSup.xlsx"
#'
#' full_path <- file.path(path_to_tox, file_name)
#' 
#' tox_list <- create_toxEval(full_path)
#' \dontrun{
#' ACClong <- get_ACC(tox_list$chem_info$CAS)
#' ACClong <- remove_flags(ACClong)
#' 
#' cleaned_ep <- clean_endPoint_info(endPointInfo)
#' filtered_ep <- filter_groups(cleaned_ep)
#' chemicalSummary <- get_chemical_summary(tox_list, ACClong, filtered_ep)
#'                                        
#' plot_tox_stacks(chemicalSummary, tox_list$chem_site, "Biological")   
#' plot_tox_stacks(chemicalSummary, tox_list$chem_site, "Chemical Class")
#' plot_tox_stacks(chemicalSummary, tox_list$chem_site, "Chemical", include_legend = FALSE) 
#' }
plot_tox_stacks <- function(chemicalSummary, 
                            chem_site,
                            category = "Biological",
                            mean_logic = FALSE,
                            manual_remove = NULL,
                            include_legend = TRUE, 
                            font_size = NA,
                            title = NA){
  
  match.arg(category, c("Biological","Chemical Class","Chemical"))
  
  site <- EAR <- sumEAR <- meanEAR <- groupCol <- nonZero <- ".dplyr"
  SiteID <- site_grouping <- n <- index <- `Short Name` <- count <- x <- y <- label <- ".dplyr"
    
  if(!("site_grouping" %in% names(chem_site))){
    chem_site$site_grouping <- ""
  }
  
  if(include_legend & category == "Chemical"){
    graphData <- graph_chem_data(chemicalSummary = chemicalSummary,
                           manual_remove = manual_remove,
                           mean_logic = mean_logic)   
    names(graphData)[names(graphData) == "maxEAR"] <- "meanEAR"
    names(graphData)[names(graphData) == "chnm"] <- "category"
  } else {
    graphData <- tox_boxplot_data(chemicalSummary = chemicalSummary,
                           category = category,
                           manual_remove = manual_remove,
                           mean_logic = mean_logic) 
    if(category == "Chemical"){
      graphData$category <- graphData$chnm
    } 
  }

  counts <- chemicalSummary %>%
    select(site, date) %>%
    distinct() %>%
    group_by(site) %>%
    summarize(count = n()) %>%
    left_join(select(chem_site, site=SiteID, `Short Name`, site_grouping), by="site") %>%
    select(-site)

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
    
    placement <- -0.05*diff(range(graphData$meanEAR))
    
    label_samples <- data.frame(x=-Inf,
                                y=placement,
                                label="# Samples", 
                                site_grouping = NA,
                                stringsAsFactors = FALSE)
    if(isTRUE(is.null(levels(chem_site$site_grouping)))){
      x <- factor(chem_site$site_grouping)
      label_samples$site_grouping <- levels(x)[1]
    } else {
      label_samples$site_grouping <- factor(levels(chem_site$site_grouping)[1],
                                            levels = levels(chem_site$site_grouping))
    }
    
    upperPlot <- ggplot(graphData, 
                        aes(x=`Short Name`, y=meanEAR, fill = category)) +
      theme_minimal() +
      xlab("") +
      ylab(paste(ifelse(mean_logic,"Mean","Maximum"), "EAR Per Site")) +
      facet_grid(. ~ site_grouping, scales="free", space="free") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
      geom_text(data = counts, 
                aes(label = count, x=`Short Name`,y = placement), 
                size=ifelse(is.na(font_size),3,0.30*font_size),inherit.aes = FALSE) +
      geom_text(data = label_samples,hjust=0.9,
                aes(x=x,y=y,label=label),
                size=ifelse(is.na(font_size),3,0.30*font_size),inherit.aes = FALSE)

  } else {

    graphData <- chemicalSummary %>%
      select(-site) 
    
    placement <- -0.05*diff(range(graphData$meanEAR))
    
    dates <- arrange(distinct(select(graphData, date))) 
    dates$index <- 1:(nrow(dates))
    
    graphData <- graphData %>%
      left_join(dates, by="date")

    if(category == "Biological"){
      graphData$category <- graphData$Bio_category
    } else if (category == "Chemical Class"){
      graphData$category <- graphData$Class
    } else {
      graphData$category <- graphData$chnm
    }
    
    upperPlot <- ggplot(graphData, aes(x=index, y=EAR, fill = category)) +
      theme_minimal() +
      theme(axis.text.x=element_blank(),
            axis.ticks.x=element_blank()) +
      xlab("Individual Samples") +
      ylab("EAR") 

  }
  
  upperPlot <- upperPlot +
    geom_col()  +
    theme(plot.margin = unit(c(5.5,5.5,5.5,12), "pt"))
  
  if(length(unique(graphData$category)) <= length(cbValues)){
    upperPlot <- upperPlot + 
      scale_fill_manual(name = category,values = cbValues, drop=TRUE)

  } 
  
  if(!include_legend){
    upperPlot <- upperPlot +
      theme(legend.position="none")
  }
  
  if(!is.na(font_size)){
    upperPlot <- upperPlot +
      theme(axis.text = element_text(size = font_size),
            strip.text = element_text(size = font_size))
  }
  
  if(!is.na(title)){
    upperPlot <- upperPlot +
      ggtitle(title)
    
    if(!is.na(font_size)){
      upperPlot <- upperPlot +
        theme(plot.title = element_text(size=font_size))
    }
  }
  
  return(upperPlot)
}