#' plot_tox_heatmap
#' 
#' The \code{plot_tox_heatmap} function creates a heat (tile) map with sites on the x-axis, 
#' a specified grouping on the y-axis (defined by the category argument), and color shading 
#' defining the mean or maximum EAR. See "Summarizing the data" in the Introduction vignette:  
#' \href{../doc/Introduction.html#summarize_data}{\code{vignette("Introduction", package = "toxEval")}} for a description on how the 
#' EAR values are computed, aggregated, and summarized. The y-axis grouping can be "Biological",
#' "Chemical Class", or "Chemical". When specifying the "Chemical" option, a secondary y-axis 
#' is automatically included to group chemicals into chemical class. The function computes 
#' default breaks for the color scale to match the spread of the data, but breaks can also 
#' be customized with the breaks argument.

#' This is a function where it may be ideal to create a custom order to the sites 
#' (for example, west-to-east). See the above section "Custom configuration"
#' \href{../doc/Introduction.html#custom_config}{\code{vignette("Introduction", package = "toxEval")}} for instructions on how to convert 
#' the character vector sites to a factor with ordered levels.
#' 
#' @param chemicalSummary Data frame from \code{\link{get_chemical_summary}}.
#' @param chem_site Data frame with columns SiteID, site_grouping, and Short Name.
#' @param category Either "Biological", "Chemical Class", or "Chemical".
#' @param manual_remove Vector of categories to remove.
#' @param mean_logic Logical.  \code{TRUE} displays the mean sample from each site,
#' \code{FALSE} displays the maximum sample from each site.
#' @param sum_logic Logical. \code{TRUE} sums the EARs in a specified grouping,
#' \code{FALSE} does not. \code{FALSE} may be better for traditional benchmarks as
#' opposed to ToxCast benchmarks.
#' @param plot_ND Logical. Whether or not to plot "Biological" groupings,
#' "Chemical Class" groupings, or "Chemical" that do not have any detections. 
#' @param font_size Numeric value to adjust the axis font size.
#' @param title Character title for plot. 
#' @param breaks Numerical vector to define data bins and legend breaks.
#' @export
#' @rdname plot_tox_heatmap
#' @import ggplot2
#' @importFrom stats median
#' @importFrom dplyr full_join filter mutate left_join right_join
#' @examples
#' path_to_tox <-  system.file("extdata", package="toxEval")
#' file_name <- "OWC_data_fromSup.xlsx"
#' full_path <- file.path(path_to_tox, file_name)
#' 
#' tox_list <- create_toxEval(full_path)
#' 
#' ACC <- get_ACC(tox_list$chem_info$CAS)
#' ACC <- remove_flags(ACC)
#' 
#' cleaned_ep <- clean_endPoint_info(endPointInfo)
#' filtered_ep <- filter_groups(cleaned_ep)
#' 
#' chemicalSummary <- get_chemical_summary(tox_list, ACC, filtered_ep)
#'                                         
#' #Order the site_groupings:
#' tox_list$chem_site$site_grouping <- factor(tox_list$chem_site$site_grouping,
#'               levels=c("Lake Superior",
#'               "Lake Michigan",
#'               "Lake Huron",
#'               "Lake Erie",
#'               "Lake Ontario"))
#' 
#' #Order sites:
#' sitesOrdered <- c("StLouis","Nemadji","WhiteWI","Bad","Montreal",
#' "PresqueIsle","Ontonagon","Sturgeon","Tahquamenon","Burns",
#' "IndianaHC","StJoseph","PawPaw","Kalamazoo","GrandMI",
#' "Milwaukee","Muskegon","WhiteMI","PereMarquette","Manitowoc",
#' "Manistee","Fox","Oconto","Peshtigo","Menominee",
#' "Indian","Cheboygan","Ford","Escanaba","Manistique",
#' "ThunderBay","AuSable","Rifle","Saginaw","BlackMI",
#' "Clinton","Rouge","HuronMI","Raisin","Maumee",
#' "Portage","Sandusky","HuronOH","Vermilion","BlackOH",
#' "Rocky","Cuyahoga","GrandOH","Cattaraugus","Tonawanda",
#' "Genesee","Oswego","BlackNY","Oswegatchie","Grass",
#' "Raquette","StRegis")
#' 
#' tox_list$chem_site$`Short Name` <- factor(tox_list$chem_site$`Short Name`,
#'               levels = sitesOrdered)
#'               
#' plot_tox_heatmap(chemicalSummary, 
#'                  tox_list$chem_site, 
#'                  category = "Biological",
#'                  manual_remove = "Undefined")
#' plot_tox_heatmap(chemicalSummary, tox_list$chem_site, category = "Chemical Class")
#' plot_tox_heatmap(chemicalSummary, tox_list$chem_site, category = "Chemical")
#' 
plot_tox_heatmap <- function(chemicalSummary, 
                             chem_site, 
                             category = "Biological",
                             breaks = c(0.00001,0.0001,0.001,0.01,0.1,1,5),
                             manual_remove = NULL,
                             mean_logic = FALSE,
                             sum_logic = TRUE,
                             plot_ND = TRUE, 
                             font_size = NA,
                             title = NA){
  
  match.arg(category, c("Biological","Chemical Class","Chemical"))
  
  SiteID <- site_grouping <- `Short Name` <- chnm <- maxEAR <- ".dplyr"
  site <- EAR <- sumEAR <- meanEAR <- ".dplyr"

  fill_label <- ifelse(mean_logic, "Mean EAR", "Max EAR")
  
  if(!("site_grouping" %in% names(chem_site))){
    chem_site$site_grouping <- "Sites"
  }
  
  if(!plot_ND){
    chemicalSummary <- chemicalSummary[chemicalSummary$EAR > 0,]
  }
  
  if(category == "Chemical"){
    plot_back <- plot_heat_chemicals(chemicalSummary=chemicalSummary, 
                                     mean_logic=mean_logic,
                                     sum_logic=sum_logic,
                                     chem_site=chem_site)
    
  } else {
    
    graphData <- tox_boxplot_data(chemicalSummary = chemicalSummary,
                           category = category,
                           manual_remove = manual_remove,
                           mean_logic = mean_logic,
                           sum_logic = sum_logic)

    graphData <- graphData %>%
      left_join(chem_site[, c("SiteID", "site_grouping", "Short Name")],
                by=c("site"="SiteID"))
    
    
    plot_back <- ggplot(data = graphData) +
      geom_tile(aes(x = `Short Name`, y=category, fill=meanEAR)) +
      theme_bw() +
      theme(axis.text.x = element_text( angle = 90,vjust=0.5,hjust = 0.975)) +
      ylab("") +
      xlab("") +
      labs(fill=fill_label) +
      scale_fill_gradient( guide = "legend",
                           trans = 'log',
                           low = "white", high = "steelblue",
                           breaks=breaks,
                           na.value = 'transparent',labels=fancyNumbers2) +
      facet_grid(. ~ site_grouping, scales="free", space="free") +
      theme(strip.text.y = element_text(angle=0, hjust=0), 
            strip.background = element_rect(fill="transparent", colour = NA),
            panel.spacing = unit(0.05, "lines"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.ticks = element_blank(),
            plot.background = element_rect(fill = "transparent",colour = NA))

  }
  
  if(!is.na(font_size)){
    plot_back <- plot_back +
      theme(axis.text = element_text(size = font_size),
            strip.text = element_text(size = font_size))
  }
  
  if(!is.na(title)){
    plot_back <- plot_back +
      ggtitle(title)
    
    if(!is.na(font_size)){
      plot_back <- plot_back +
        theme(plot.title = element_text(size=font_size))
    }
  }  
  
  return(plot_back)
}


plot_heat_chemicals <- function(chemicalSummary,
                                chem_site,
                                mean_logic,
                                sum_logic){
  
  SiteID <- site_grouping <- `Short Name` <- chnm <- maxEAR <- ".dplyr"
  site <- EAR <- sumEAR <- meanEAR <- ".dplyr"
  
  graphData <- graph_chem_data(chemicalSummary, 
                               mean_logic=mean_logic,
                               sum_logic = sum_logic)
  
  if(!("site_grouping" %in% names(chem_site))){
    chem_site$site_grouping <- "Sites"
  }
  
  graphData <- graphData %>%
    left_join(chem_site[, c("SiteID", "site_grouping", "Short Name")],
              by=c("site"="SiteID"))
  
  fill_text <- ifelse(mean_logic, "Mean EAR", "Max EAR")
  
  heat <- ggplot(data = graphData) +
    geom_tile(aes(x = `Short Name`, y=chnm, fill=meanEAR)) +
    theme_bw() +
    theme(axis.text.x = element_text( angle = 90,vjust=0.5,hjust = 1)) +
    ylab("") +
    xlab("") +
    labs(fill=fill_text) +
    scale_fill_gradient( guide = "legend",
                         trans = 'log',
                         low = "white", high = "steelblue",
                         breaks=c(0.00001,0.0001,0.001,0.01,0.1,1,5),
                         na.value = 'transparent',labels=fancyNumbers2) +
    facet_grid(Class ~ site_grouping, scales="free", space="free") +
    theme(strip.text.y = element_text(angle=0, hjust=0), 
          strip.background = element_rect(fill="transparent", colour = NA),
          # axis.text.y = element_text(face=ifelse(levels(graphData$category) %in% c("Total"),"bold","italic")),
          panel.spacing = unit(0.05, "lines"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.background = element_rect(fill = "transparent",colour = NA))
  
  return(heat)
  
}



