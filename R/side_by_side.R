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
#' 
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