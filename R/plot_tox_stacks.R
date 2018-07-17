#' Plot stacked bar charts
#' 
#' The \code{plot_tox_stacks} function creates a set of boxplots representing EAR 
#' values computed with the \code{\link{get_chemical_summary}} function, and 
#' dependent on the choice of several input options. See "Summarizing the data" 
#' in the Introduction vignette: \href{../doc/Introduction.html#summarize_data}{\code{vignette("Introduction", package = "toxEval")}}
#' for a description on how the EAR values are computed, aggregated, and summarized. 
#' Choosing "Chemical Class" in the category argument will generate separate stacked 
#' bars for each unique class. "Chemical" will generate stacked bars for each individual 
#' chemical, and "Biological" will generate stacked bars for each group in the selected 
#' ToxCast annotation. The legend can optionally be turned on or off using the
#' include_legend argument. It may be impractical for instance to show the 
#' legend for "Chemical" if there are hundreds of chemicals.
#' 
#' The graph displays a slightly different result for a single site. Providing 
#' data with only one site displays each individual sample as a stacked bar 
#' rather than the mean or maximum for a site. 
#' 
#' This is a function where it may be ideal to create a custom order to the sites 
#' (for example, west-to-east). See the above section "Custom configuration"
#' \href{../doc/Introduction.html#custom_config}{\code{vignette("Introduction", package = "toxEval")}} for instructions on how to convert 
#' the character vector sites to a factor with ordered levels.
#' 
#' @param chemicalSummary Data frame from \code{\link{get_chemical_summary}}.
#' @param category Character. Either "Biological", "Chemical Class", or "Chemical".
#' @param chem_site Data frame with at least columns SiteID, site_grouping, and Short Name.
#' @param mean_logic Logical.  \code{TRUE} displays the mean sample from each site,
#' \code{FALSE} displays the maximum sample from each site.
#' @param sum_logic Logical. \code{TRUE} sums the EARs in a specified grouping,
#' \code{FALSE} does not. \code{FALSE} may be better for traditional benchmarks as
#' opposed to ToxCast benchmarks.
#' @param manual_remove Vector of categories to remove.
#' @param include_legend Logical. Used to include legend or not.
#' @param font_size Numeric value to adjust the axis font size.
#' @param title Character title for plot. 
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
#' 
#' ACC <- get_ACC(tox_list$chem_info$CAS)
#' ACC <- remove_flags(ACC)
#' 
#' cleaned_ep <- clean_endPoint_info(end_point_info)
#' filtered_ep <- filter_groups(cleaned_ep)
#' chemicalSummary <- get_chemical_summary(tox_list, ACC, filtered_ep)
#'                                        
#' plot_tox_stacks(chemicalSummary, tox_list$chem_site, "Biological")   
#' plot_tox_stacks(chemicalSummary, tox_list$chem_site, "Chemical Class")
#' plot_tox_stacks(chemicalSummary, tox_list$chem_site, "Chemical", include_legend = FALSE) 
#' 
plot_tox_stacks <- function(chemicalSummary, 
                            chem_site,
                            category = "Biological",
                            mean_logic = FALSE,
                            sum_logic = TRUE,
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
  
  if(category == "Chemical"){
    graphData <- graph_chem_data(chemicalSummary = chemicalSummary,
                           manual_remove = manual_remove,
                           mean_logic = mean_logic,
                           sum_logic = sum_logic)   
    names(graphData)[names(graphData) == "maxEAR"] <- "meanEAR"
    names(graphData)[names(graphData) == "chnm"] <- "category"
  } else {
    graphData <- tox_boxplot_data(chemicalSummary = chemicalSummary,
                           category = category,
                           manual_remove = manual_remove,
                           mean_logic = mean_logic,
                           sum_logic = sum_logic) 
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
  single_site <- length(siteToFind) == 1
  
  if(!single_site){
    
    y_label <- fancyLabels(category, mean_logic, sum_logic, single_site, sep = TRUE)
    
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
      ylab(y_label[["y_label"]]) +
      facet_grid(. ~ site_grouping, scales="free", space="free") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
      geom_text(data = counts, 
                aes(label = count, x=`Short Name`,y = placement), 
                size=ifelse(is.na(font_size),3,0.30*font_size),inherit.aes = FALSE) +
      geom_text(data = label_samples,hjust=1,
                aes(x=x,y=y,label=label),
                size=ifelse(is.na(font_size),2,0.25*font_size),inherit.aes = FALSE) +
      labs(caption = y_label[["caption"]])  

  } else {

    y_label <- "EARs per Individual Sample"
    
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
      ylab(y_label) 
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
            strip.text = element_text(size = font_size),
            axis.title =   element_text(size=font_size))
  }
  
  if(!is.na(title)){
    upperPlot <- upperPlot +
      ggtitle(title)
    
    if(!is.na(font_size)){
      upperPlot <- upperPlot +
        theme(plot.title = element_text(size=font_size),
              plot.caption = element_text(size=font_size))
    }
  }
  
  if(packageVersion("ggplot2") >= '3.0.0'){
    upperPlot <- upperPlot +
      coord_cartesian(clip = "off")
  } 

  return(upperPlot)
}

