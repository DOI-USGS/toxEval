#' Summary of hits per site/category
#' 
#' These functions create a table with several rows per site depending on which 
#' categories get hits. A "hit" is based on a user specified hit_threshold. For example, 
#' if "Biological" is the category, and a site has hits above a threshold for 
#' "DNA Binding" and "Nuclear Receptors", that site will have 2 rows of data in this table.
#'
#' For each row, there are 4 columns. Site and category (as defined by the category argument) 
#' define the row. "Hits per Sample" are how many samples at that site had hits (based on 
#' hit_threshold) for the category. "Number of Samples" is number of samples.
#' 
#' The tables show slightly different results for a single site. Instead of one row per 
#' site/category, there is one row per category.
#' 
#' @param chemicalSummary data frame from \code{get_chemical_summary}
#' @param mean_logic logical \code{TRUE} is mean, \code{FALSE} is maximum
#' @param category either "Biological", "Chemical Class", or "Chemical"
#' @param hit_threshold numeric threshold defining a "hit"
#' @export
#' @rdname hits_summary_DT
#' @import DT
#' @importFrom stats median
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
#' stats_group <- hits_summary(chemicalSummary, "Biological")
#'
#' hits_summary_DT(chemicalSummary, category = "Biological")
#' hits_summary_DT(chemicalSummary, category = "Chemical Class")
#' hits_summary_DT(chemicalSummary, category = "Chemical")
#' }
hits_summary_DT <- function(chemicalSummary, 
                          category = "Biological",
                          mean_logic = FALSE,
                          hit_threshold = 0.1){
    
  chnm <- Class <- Bio_category <- site <- EAR <- sumEAR <- hits <- n <- ".dplyr"
  
  match.arg(category, c("Biological","Chemical Class","Chemical"))
  
  hits_summaryOrdered <- hits_summary(chemicalSummary = chemicalSummary,
               category = category,
               hit_threshold = hit_threshold, 
               mean_logic=mean_logic)
  
  siteToFind <- unique(chemicalSummary$site)

  meanChem <- grep("Samples with hits",names(hits_summaryOrdered))
  
  colToSort <- meanChem-1

  tableGroup <- DT::datatable(hits_summaryOrdered, extensions = 'Buttons',
                              rownames = FALSE,
                              options = list(order=list(list(colToSort,'desc')),
                                            dom = 'Bfrtip',
                                            buttons = list('colvis')))
  
  tableGroup <- formatStyle(tableGroup, names(hits_summaryOrdered)[meanChem],
                            background = styleColorBar(range(hits_summaryOrdered[,meanChem],na.rm=TRUE), 'wheat'),
                            backgroundSize = '100% 90%',
                            backgroundRepeat = 'no-repeat',
                            backgroundPosition = 'center')

  return(tableGroup)
}

#' @export
#' @rdname hits_summary_DT
#' @return data frame with columns "Hits per Sample", "Individual Hits",
#' "nSample", "site", and "category"
hits_summary <- function(chemicalSummary, category, 
                         hit_threshold = 0.1, mean_logic = FALSE){
  
  Class <- Bio_category <- `Hits per Sample` <- site <- EAR <- sumEAR <- chnm <- n <- hits <- ".dplyr"
  
  siteToFind <- unique(chemicalSummary$site)
  
  if(category == "Chemical"){
    chemicalSummary <- mutate(chemicalSummary, category = chnm)
  } else if (category == "Chemical Class"){
    chemicalSummary <- mutate(chemicalSummary, category = Class)
  } else {
    chemicalSummary <- mutate(chemicalSummary, category = Bio_category)
  }
  
  chemicalSummary <- select(chemicalSummary, -Class, -Bio_category, -chnm)
  
  if(length(siteToFind) == 1){
    chemicalSummary$site <- chemicalSummary$category
  } else {
    chemicalSummary$site <- chemicalSummary$shortName
  }
  
  hits_summary <- chemicalSummary %>%
    group_by(site, date,category) %>%
    summarise(sumEAR = sum(EAR),
              hits = sum(EAR > hit_threshold)) %>%
    group_by(site,category) %>%
    summarise(`Samples with hits` = sum(sumEAR > hit_threshold),
              nSamples = n()) %>%
    arrange(desc(`Samples with hits`))
  
  if(length(siteToFind) == 1){
    hits_summary <- hits_summary[,c("category","Samples with hits","nSamples")]
  }
  
  hits_summary <- rename(hits_summary, `Number of Samples`=nSamples)
  
  return(hits_summary)
}