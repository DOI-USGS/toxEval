#' Plot EAR heat maps
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
#' If there are site/parameters (chemical/chemical class/biological grouping) combinations that
#' don't have data, those areas are represented by an "X". If there are 0 values,
#' they are considered "non-detects", and represented with a distinct color.
#' 
#' @param chemical_summary Data frame from \code{\link{get_chemical_summary}}.
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
#' @examples
#' path_to_tox <- system.file("extdata", package="toxEval")
#' file_name <- "OWC_data_fromSup.xlsx"
#' full_path <- file.path(path_to_tox, file_name)
#' 
#' tox_list <- create_toxEval(full_path)
#' 
#' ACC <- get_ACC(tox_list$chem_info$CAS)
#' ACC <- remove_flags(ACC)
#' 
#' cleaned_ep <- clean_endPoint_info(end_point_info)
#' filtered_ep <- filter_groups(cleaned_ep)
#' 
#' chemical_summary <- get_chemical_summary(tox_list, ACC, filtered_ep)
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
#' plot_tox_heatmap(chemical_summary, 
#'                  tox_list$chem_site, 
#'                  category = "Biological",
#'                  manual_remove = "Undefined")
#' plot_tox_heatmap(chemical_summary, tox_list$chem_site, category = "Chemical Class")
#' plot_tox_heatmap(chemical_summary, tox_list$chem_site, category = "Chemical")
#' 
plot_tox_heatmap <- function(chemical_summary, 
                             chem_site, 
                             category = "Biological",
                             breaks = c(0.00001,0.0001,0.001,0.01,0.1,1,10),
                             manual_remove = NULL,
                             mean_logic = FALSE,
                             sum_logic = TRUE,
                             plot_ND = TRUE, 
                             font_size = NA,
                             title = NA){
  
  match.arg(category, c("Biological","Chemical Class","Chemical"))
  
  SiteID <- site_grouping <- `Short Name` <- chnm <- maxEAR <- ".dplyr"
  site <- EAR <- sumEAR <- meanEAR <- ".dplyr"

  if(!("site_grouping" %in% names(chem_site))){
    chem_site$site_grouping <- "Sites"
  }
  
  if(!plot_ND){
    chemical_summary <- chemical_summary[chemical_summary$EAR > 0,]
  }
  
  if(category == "Chemical"){
    plot_back <- plot_heat_chemicals(chemical_summary=chemical_summary, 
                                     mean_logic=mean_logic,
                                     sum_logic=sum_logic,
                                     chem_site=chem_site,
                                     breaks=breaks)
    
  } else {
    
    graphData <- tox_boxplot_data(chemical_summary = chemical_summary,
                           category = category,
                           manual_remove = manual_remove,
                           mean_logic = mean_logic,
                           sum_logic = sum_logic)

    graphData <- graphData %>%
      left_join(chem_site[, c("SiteID", "site_grouping", "Short Name")],
                by=c("site"="SiteID"))
    
    # This requires non-detects to be 0. If that changes we'll need to update:
    graphData$meanEAR[graphData$meanEAR == 0] <- NA
    
    complete_data_filled <- get_complete_set_category(chemical_summary, graphData, chem_site, category)
    
    any_missing <- nrow(complete_data_filled) > nrow(graphData)
    any_non_detects <- any(is.na(graphData$meanEAR))
    
    single_site <- length(unique(chemical_summary$site)) == 1
    
    y_label <- fancyLabels(category, mean_logic, sum_logic, single_site, sep = TRUE, include_site = FALSE)
    
    caption <- y_label[['caption']]
    fill_label <- y_label[['y_label']]
    
    plot_back <- ggplot(data = graphData) +
      geom_point(data = complete_data_filled, aes(x = `Short Name`, y=category, shape=""), size = 2 ) +
      geom_tile(aes(x = `Short Name`, y=category, fill=meanEAR, color="")) +
      theme_bw() +
      theme(axis.text.x = element_text( angle = 90,vjust=0.5,hjust = 0.975)) +
      ylab("") +
      xlab("") +
      scale_colour_manual(values=NA) + 
      scale_shape_manual(values=4) +
      labs(fill = fill_label, caption = caption) +
      # labs(fill=y_label[["y_label"]], caption = y_label[["caption"]]) +
      scale_fill_gradient( trans = 'log',
                           low = "white", high = "steelblue",
                           breaks=breaks,
                           na.value = 'khaki',labels=fancyNumbers2) +
      facet_grid(. ~ site_grouping, scales="free", space="free") +
      theme(strip.text.y = element_text(angle=0, hjust=0), 
            strip.background = element_rect(fill="transparent", colour = NA),
            panel.spacing = unit(0.05, "lines"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.ticks = element_blank(),
            plot.background = element_rect(fill = "transparent",colour = NA))

    any_non_detects <- any(is.na(graphData$meanEAR))
    complete_data_filled <- get_complete_set_category(chemical_summary, graphData, chem_site)
    
    if(any_non_detects & any_missing){
      plot_back <- plot_back +
        guides(colour=guide_legend("Non-detects", override.aes=list(colour="khaki", fill="khaki"), order = 2),
               shape=guide_legend("Missing", order = 3),
               fill = guide_colorbar(order=1)) 
    } else if (any_non_detects){
      plot_back <- plot_back +
        guides(colour=guide_legend("Non-detects", override.aes=list(colour="khaki", fill="khaki"), order = 2),
               fill = guide_colorbar(order=1),
               shape = "none")    
    } else if (any_missing){
      plot_back <- plot_back +
        guides(shape=guide_legend("Missing", order = 2),
               fill = guide_colorbar(order=1),
               colour = "none")
    } else {
      plot_back <- plot_back +
        guides(fill = guide_colorbar(order=1),
               shape = "none",
               colour = "none")
    }
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


plot_heat_chemicals <- function(chemical_summary,
                                chem_site,
                                mean_logic,
                                sum_logic,
                                breaks){
  
  SiteID <- site_grouping <- `Short Name` <- chnm <- maxEAR <- ".dplyr"
  site <- EAR <- sumEAR <- meanEAR <- ".dplyr"
  
  graphData <- graph_chem_data(chemical_summary, 
                               mean_logic=mean_logic,
                               sum_logic = sum_logic)
  
  if(!("site_grouping" %in% names(chem_site))){
    chem_site$site_grouping <- "Sites"
  }
  single_site <- length(unique(chemical_summary$site)) == 1

  y_label <- fancyLabels(category = "Chemical", mean_logic, sum_logic, single_site, sep = TRUE, include_site = FALSE)
  
  caption <- y_label[['caption']]
  fill_text <- y_label[['y_label']]

  graphData <- graphData %>%
    left_join(chem_site[, c("SiteID", "site_grouping", "Short Name")],
              by=c("site"="SiteID"))
  
  # This requires non-detects to be 0. If that changes we'll need to update:
  graphData$meanEAR[graphData$meanEAR == 0] <- NA
  
  complete_data_filled <- get_complete_set(chemical_summary, graphData, chem_site)
  
  any_missing <- nrow(complete_data_filled) > nrow(graphData)
  any_non_detects <- any(is.na(graphData$meanEAR))
  
  heat <- ggplot(data = graphData) +
    geom_point(data = complete_data_filled, aes(x = `Short Name`, y=chnm, shape=""), size = 2 ) +
    geom_tile(aes(x = `Short Name`, y=chnm, fill=meanEAR, color = "")) +
    theme_bw() +
    theme(axis.text.x = element_text( angle = 90,vjust=0.5,hjust = 1)) +
    ylab("") +
    xlab("") +
    labs(fill = fill_text, caption = caption) +
    scale_fill_gradient( na.value = 'khaki',
                         trans = 'log', low = "white", high = "steelblue",
                         breaks=breaks,
                         labels=fancyNumbers2, guide = "colourbar") +
    scale_colour_manual(values=NA) + 
    scale_shape_manual(values=4) +
    facet_grid(Class ~ site_grouping, scales="free", space="free") +
    theme(strip.text.y = element_text(angle=0, hjust=0), 
          strip.background = element_rect(fill="transparent", colour = NA),
          # axis.text.y = element_text(face=ifelse(levels(graphData$category) %in% c("Total"),"bold","italic")),
          panel.spacing = unit(0.05, "lines"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.background = element_rect(fill = "transparent",colour = "transparent")) 
  
  if(any_non_detects & any_missing){
    heat <- heat +
            guides(colour=guide_legend("Non-detects", override.aes=list(colour="khaki", fill="khaki"), order = 2),
                   shape=guide_legend("Missing", order = 3),
                   fill = guide_colorbar(order=1)) 
  } else if (any_non_detects){
    heat <- heat +
      guides(colour=guide_legend("Non-detects", override.aes=list(colour="khaki", fill="khaki"), order = 2),
             fill = guide_colorbar(order=1),
             shape = "none")    
  } else if (any_missing){
    heat <- heat +
      guides(shape=guide_legend("Missing", order = 2),
             fill = guide_colorbar(order=1),
             colour = "none")
  } else {
    heat <- heat +
      guides(fill = guide_colorbar(order=1),
             shape = "none",
             colour = "none")
  }
  
  return(heat)
  
}

# There's probably a faster way to do this:
get_complete_set <- function(chemical_summary, graphData, chem_site){
  
  `Short Name` <- site_grouping <- Class <- chnm <- ".dplyr"
  
  complete_data <- select(chem_site, `Short Name`, site_grouping)
  complete_data_filled <- data.frame()
  
  for(chms in levels(chemical_summary$chnm)){
    complete_data$chnm <- chms
    complete_data_filled <- bind_rows(complete_data_filled, complete_data)
  }
  complete_data_filled$chnm <- factor(complete_data_filled$chnm, levels = levels(graphData$chnm))
  complete_data_filled <- left_join(complete_data_filled, distinct(select(chemical_summary, chnm, Class)), by="chnm")
  return(complete_data_filled)
}

# There's probably a faster way to do this:
get_complete_set_category <- function(chemical_summary, graphData, chem_site, category){
  
  `Short Name` <- site_grouping <- Class <- chnm <- ".dplyr"
  
  complete_data <- select(chem_site, `Short Name`, site_grouping)
  complete_data_filled <- data.frame()
  categories <- levels(graphData$category)
  
  for(cats in categories){
    complete_data$category <- cats
    complete_data_filled <- bind_rows(complete_data_filled, complete_data)
  }
  complete_data_filled$category <- factor(complete_data_filled$category, levels = levels(graphData$category))

  return(complete_data_filled)
}