#' Prepare boxplot data
#' 
#' A set of functions to prepare the data for boxplots. Often, these
#' functions are used within the plotting functions. They are exported however
#' to allow custom graphs to be created. 
#' 
#' The function side_by_side_data will combine two data frames,
#' either the output of \code{\link{get_chemical_summary}} or \code{\link{graph_chem_data}},
#' into a single data frame. The important work here is that the chemicals
#' and classes factor levels are ordered primarily based on "gd_left", but
#' include "gd_right" when the contents are mismatched.
#' 
#' @export
#' @rdname graph_data_prep
#' 
#' @param gd_left Data frame that must include the columns chnm, Class, and either EAR or meanEAR.
#' @param gd_right Data frame that must include the columns chnm, Class, and either EAR or meanEAR.
#' @param left_title Character that will be associated with the "gd_left" data 
#' frame in a column named "guide_side". 
#' @param right_title Character that will be associated with the "gd_right" data 
#' frame in a column named "guide_side". 
#' 
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
#' cleaned_ep <- clean_endPoint_info(end_point_info)
#' filtered_ep <- filter_groups(cleaned_ep)
#' 
#' chemical_summary <- get_chemical_summary(tox_list, ACC, filtered_ep)
#' # Let's say we want to compare 2 chemical summaries
#' # We'll look at one summing EARs, and with concentrations
#' # First, we need a chemical summary for concentrations:
#' chemical_summary_conc <- get_concentration_summary(tox_list)
#' 
#' gd_tox <- graph_chem_data(chemical_summary)
#' gd_conc <- graph_chem_data(chemical_summary_conc)
#' 
#' ch_combo <- side_by_side_data(gd_tox, gd_conc, 
#'                               left_title = "ToxCast", 
#'                               right_title = "Concentrations")
#' plot_chemical_boxplots(ch_combo, guide_side, 
#'                        x_label = "") +
#'   ggplot2::facet_grid(. ~ guide_side, scales = "free_x")
side_by_side_data <- function(gd_left, 
                              gd_right,
                              left_title = "Left",
                              right_title = "Right"){
  
  Class <- EAR <- chnm <- hit_label <- hits <- meanEAR <- ymax <- ymin <- ".dplyr"
  
  gd_left$guide_side <- left_title
  gd_right$guide_side <- right_title
  
  gd_left_no_factor <- gd_left %>% 
    mutate(chnm = as.character(chnm),
           Class = as.character(Class))
  
  gd_right_no_factor <- gd_right %>% 
    mutate(chnm = as.character(chnm),
           Class = as.character(Class))
  
  chem_data_no_factors <- gd_left_no_factor %>% 
    bind_rows(gd_right_no_factor %>% 
                filter(!(chnm %in% unique(gd_left_no_factor$chnm))))
  
  graph_data <- "meanEAR" %in% names(chem_data_no_factors)
  
  if(!(graph_data)){
    chem_data_no_factors <- chem_data_no_factors %>% 
      rename(meanEAR = EAR)
  }
  
  orderChem_1_2 <- chem_data_no_factors %>%
    group_by(chnm, Class) %>%
    summarise(median = quantile(meanEAR[meanEAR != 0],0.5)) %>%
    ungroup()
  
  class_order <- orderClass(chem_data_no_factors)
  
  orderChem_1_2 <- orderChem_1_2 %>%
    mutate(Class = factor(Class, levels=class_order$Class)) %>%
    arrange(desc(Class), desc(!is.na(median)), median)
  
  gd_1_2 <- bind_rows(gd_left_no_factor, 
                      gd_right_no_factor)
  gd_1_2$Class <- factor(gd_1_2$Class, levels = class_order$Class)
  gd_1_2$chnm <- factor(gd_1_2$chnm, levels = orderChem_1_2$chnm)
  
  gd_1_2$guide_side <- factor(gd_1_2$guide_side, 
                              levels = c(left_title, right_title))
  
  return(gd_1_2)
  
}

#' Create concentration summary
#' 
#' Use this function to create a chemical_summary, but instead
#' of using any benchmarks, the EAR column is simply
#' the concentration. The output of this function can be used
#' in any of the plotting or table functions in the same way 
#' that the output of \code{\link{get_chemical_summary}}.
#' 
#' 
#' @param tox_list List with data frames for chem_data, chem_info, and chem_site.
#' Created with \code{\link{create_toxEval}}.
#' @param chem_data \emph{Optional} data frame with (at least) columns: CAS, SiteID, and Value. Default is \code{NULL}. 
#' The argument will over-ride what is in tox_list.
#' @param chem_site \emph{Optional} data frame with (at least) columns: SiteID, and Short Name. Default is \code{NULL}. 
#' The argument will over-ride what is in tox_list.
#' @param chem_info \emph{Optional} data frame with (at least) columns: CAS, and class. Default is \code{NULL}. 
#' The argument will over-ride what is in tox_list.
#' @param tox_names Logical whether to use the provided chemical names from the ToxCast or not. If 
#' there is not a match by CAS, the function will look for a column "Chemical" in the "Chemical"
#' tab. If that column doesn't exist, it will create a (not good!) name.
#' @export
#' @return a data frame with the columns: CAS, chnm (chemical name
#' as a factor), site, date, EAR (which is just concentration), Bio_category, shortName (of site), Class. The output of this 
#' function is where you find EAR values for every chemical/endpoint combination.
#' 
#' @examples
#' path_to_tox <-  system.file("extdata", package="toxEval")
#' file_name <- "OWC_data_fromSup.xlsx"
#' full_path <- file.path(path_to_tox, file_name)
#' 
#' tox_list <- create_toxEval(full_path)
#' 
#' chemical_summary_conc <- get_concentration_summary(tox_list)
#' head(chemical_summary_conc)  
#' plot_tox_boxplots(chemical_summary_conc, 
#'                   category = "Chemical",
#'                   x_label = "Concentration [ug/L]")                        
get_concentration_summary <- function(tox_list, 
                                      chem_data=NULL, 
                                      chem_site=NULL, 
                                      chem_info=NULL,
                                      tox_names = TRUE){
  # Getting rid of NSE warnings:
  chnm <- endPoint <- ACC_value <- Substance_Name <- Value <- `Sample Date` <- SiteID <- ".dplyr"
  EAR <- `Short Name` <- Substance_CASRN <- CAS <- Chemical <- Class <- site <- casrn <- groupCol <- ".dplyr"
  
  if(is.null(chem_data)){
    chem_data <- tox_list[["chem_data"]]
  } else {
    chem_data <- rm_em_dash(chem_data)
  }
  
  if(is.null(chem_site)){
    chem_site <- tox_list[["chem_site"]]
  } else {
    chem_site <- rm_em_dash(chem_site)
  }
  
  if(is.null(chem_info)){
    chem_info <- tox_list[["chem_info"]]
  } else {
    chem_info <- rm_em_dash(chem_info)
  }
  
  if(is.character(chem_data$Value)){
    chem_data$Value <- as.numeric(chem_data$Value)
  }
  
  chemical_summary <- chem_data %>% 
    select(CAS, SiteID, Value, `Sample Date`) %>%
    filter(!is.na(Value)) %>%
    rename(EAR = Value,
           site = SiteID,
           date = `Sample Date`) %>%
    mutate(Bio_category = "Concentration",
           endPoint = "Concentration") %>% 
    left_join(select(chem_site, 
                     site=SiteID, 
                     shortName = `Short Name`),
              by="site") %>%
    left_join(select(chem_info, CAS, Class), by="CAS") 
  
  if(tox_names){
    tox_names_key <- tox_chemicals %>% 
      select(CAS = Substance_CASRN, chnm = Substance_Name) %>% 
      filter(CAS %in% chem_info$CAS)
    
    chemical_summary <- chemical_summary %>% 
      left_join(tox_names_key, by = "CAS")
    
    if(any(is.na(chemical_summary$chnm))){
      if("Chemical" %in% names(chem_info)){
        chemical_summary <- chemical_summary %>% 
          left_join(select(chem_info, CAS, Chemical), by = "CAS")        
      } else {
        chemical_summary$Chemical = as.character(as.numeric(factor(chemical_summary$CAS)))
      }
      
      chemical_summary$chnm[is.na(chemical_summary$chnm)] <- chemical_summary$Chemical[is.na(chemical_summary$chnm)]
      
      chemical_summary <- select(chemical_summary, -Chemical)
      
    }
    
  } else {
    if("Chemical" %in% names(chem_info)){
      chemical_summary <- chemical_summary %>% 
        left_join(select(chem_info, CAS, chnm = Chemical), by = "CAS")        
    } else {
      message("Add a Chemical column to the Chemicals tab to get custom chemical names")
      chemical_summary$chnm = as.character(as.numeric(factor(chemical_summary$CAS)))
    }    
  }
    
  graphData <- graph_chem_data(chemical_summary)
  
  orderClass_df <- orderClass(graphData)
  
  orderChem_df <- orderChem(graphData, orderClass_df)
  
  chemical_summary$chnm <- factor(chemical_summary$chnm,
                                  levels = unique(orderChem_df$chnm))
  
  chemical_summary$Class <- factor(chemical_summary$Class,
                                   levels = rev(levels(orderChem_df$Class)))
  
  return(chemical_summary)
  

}
