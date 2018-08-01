#' Summary of hits per site/category
#' 
#' The \code{hits_summary_DT} (DT option) and \code{hits_summary} (data frame 
#' option) functions create tables information on the number of 
#' \code{hit_threshold} exceedances per site for each individual grouping. The 
#' table has one row per group per site that has \code{hit_threshold} 
#' exceedances. For example, if "Biological" is the category, and a 
#' site has EAR levels above the specified \code{hit_threshold} for "DNA 
#' Binding" and "Nuclear Receptors", that site will have 2 rows of data in this 
#' table.
#' 
#' For each row, there are 4 columns. Site and category (as defined by the category 
#' argument) define the row. "Samples with hits" are how many samples exceeded the 
#' hit_threshold for the specified category at the specified site. "Number of 
#' Samples" indicates how many samples were collected at an individual site 
#' based on unique date.
#' 
#' The tables contain slightly different results for evaluation of a single site. 
#' There are three columns (the Site column is dropped), and rather than one row 
#' per site/category, there is one row per category.
#' 
#' @param chemical_summary Data frame from \code{\link{get_chemical_summary}}.
#' @param sum_logic Logical. \code{TRUE} sums the EARs in a specified grouping,
#' \code{FALSE} does not. \code{FALSE} may be better for traditional benchmarks as
#' opposed to ToxCast benchmarks.
#' @param category Character. Either "Biological", "Chemical Class", or "Chemical".
#' @param hit_threshold Numeric threshold defining a "hit".
#' @export
#' @rdname hits_summary_DT
#' @return data frame with with one row per unique site/category combination. The columns
#' are site, category, Samples with Hits, and Number of Samples.
#' @importFrom stats median
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
#' chemical_summary <- get_chemical_summary(tox_list, ACC, filtered_ep)
#'
#' stats_group <- hits_summary(chemical_summary, "Biological")
#' 
#' hits_summary_DT(chemical_summary, category = "Biological")
#' hits_summary_DT(chemical_summary, category = "Chemical Class")
#' hits_summary_DT(chemical_summary, category = "Chemical")
#' 
hits_summary_DT <- function(chemical_summary, 
                          category = "Biological",
                          sum_logic = TRUE,
                          hit_threshold = 0.1){
    
  chnm <- Class <- Bio_category <- site <- EAR <- sumEAR <- hits <- n <- `Samples with
  hits` <- ".dplyr"
  
  match.arg(category, c("Biological","Chemical Class","Chemical"))

  hits_summaryOrdered <- hits_summary(chemical_summary = chemical_summary,
               category = category,
               hit_threshold = hit_threshold, 
               sum_logic = sum_logic)
  
  siteToFind <- unique(chemical_summary$site)

  meanChem <- grep("Samples with hits",names(hits_summaryOrdered))
  
  colToSort <- meanChem-1

  tableGroup <- DT::datatable(hits_summaryOrdered, extensions = 'Buttons',
                              rownames = FALSE,
                              options = list(order=list(list(colToSort,'desc')),
                                            dom = 'Bfrtip',
                                            buttons = list('colvis')))
  
  tableGroup <- DT::formatStyle(tableGroup, names(hits_summaryOrdered)[meanChem],
                            background = DT::styleColorBar(range(hits_summaryOrdered[,meanChem],na.rm=TRUE), 'wheat'),
                            backgroundSize = '100% 90%',
                            backgroundRepeat = 'no-repeat',
                            backgroundPosition = 'center')

  return(tableGroup)
}

#' @export
#' @rdname hits_summary_DT
#' @return data frame with columns "Hits per Sample", "Individual Hits",
#' "nSample", "site", and "category"
hits_summary <- function(chemical_summary, 
                         category, 
                         hit_threshold = 0.1, 
                         sum_logic = TRUE){
  
  Class <- Bio_category <- `Samples with hits` <- nSamples <- site <- EAR <- sumEAR <- chnm <- n <- hits <- ".dplyr"

  siteToFind <- unique(chemical_summary$site)

  if(category == "Chemical"){
    chemical_summary <- dplyr::mutate(chemical_summary, category = chnm)
  } else if (category == "Chemical Class"){
    chemical_summary <- dplyr::mutate(chemical_summary, category = Class)
  } else {
    chemical_summary <- dplyr::mutate(chemical_summary, category = Bio_category)
  }
  
  chemical_summary <- dplyr::select(chemical_summary, -Class, -Bio_category, -chnm)
  
  if(length(siteToFind) == 1){
    chemical_summary$site <- chemical_summary$category
  } else {
    chemical_summary$site <- chemical_summary$shortName
  }
  
  if(!sum_logic){
    hits_summary <- chemical_summary %>%
      dplyr::group_by(site,category,date) %>%
      dplyr::summarise(hits = sum(EAR > hit_threshold)) %>%
      dplyr::group_by(site,category) %>%
      dplyr::summarise(`Samples with hits` = sum(hits >= 1),
                nSamples = n()) %>%   
      dplyr::arrange(dplyr::desc(`Samples with hits`))    
  } else {
    hits_summary <- chemical_summary %>%
      dplyr::group_by(site, date,category) %>%
      dplyr::summarise(sumEAR = sum(EAR),
                hits = sum(EAR > hit_threshold)) %>%
      dplyr::group_by(site,category) %>%
      dplyr::summarise(`Samples with hits` = sum(sumEAR > hit_threshold),
                nSamples = n()) %>%
      dplyr::arrange(dplyr::desc(`Samples with hits`))    
  }

  
  if(length(siteToFind) == 1){
    hits_summary <- hits_summary[,c("category","Samples with hits","nSamples")]
  }
  
  hits_summary <- dplyr::rename(hits_summary, `Number of Samples`=nSamples)
  
  return(hits_summary)
}