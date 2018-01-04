#' table_tox_sum
#' 
#' Table of sums
#' @param chemicalSummary data frame from \code{get_chemical_summary}
#' @param mean_logic logical \code{TRUE} is mean, \code{FALSE} is maximum
#' @param category either "Biological", "Chemical Class", or "Chemical"
#' @param hit_threshold numeric threshold defining a "hit"
#' @export
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
#' chemicalSummary <- get_chemical_summary(ACClong, filtered_ep,
#'                                         tox_list)
#' }
#' # The example workflow takes a bit of time to load and compute, 
#' # so an example chemicalSummary is included pre-calculated in the package. 
#' chemicalSummary <- ex_chemSum #loading example data
#' table_tox_sum(chemicalSummary, category = "Biological")
#' table_tox_sum(chemicalSummary, category = "Chemical Class")
#' table_tox_sum(chemicalSummary, category = "Chemical")
table_tox_sum <- function(chemicalSummary, 
                          category = "Biological",
                          mean_logic = FALSE,
                          hit_threshold = 0.1){
    
  chnm <- Class <- Bio_category <- site <- EAR <- sumEAR <- hits <- n <- ".dplyr"
  
  match.arg(category, c("Biological","Chemical Class","Chemical"))
  
  statsOfGroupOrdered <- statsOfGroup(chemicalSummary = chemicalSummary,
               category = category,
               hit_threshold = hit_threshold)
  siteToFind <- unique(chemicalSummary$site)

  meanChem <- grep("Individual.Hits",names(statsOfGroupOrdered))
  maxChem <- grep("Hits.per.Sample",names(statsOfGroupOrdered))
  
  colToSort <- maxChem-1

  tableGroup <- DT::datatable(statsOfGroupOrdered,  extensions = 'Buttons',
                              rownames = FALSE,
                              options = list(pageLength = nrow(statsOfGroupOrdered),
                                            order=list(list(colToSort,'desc')),
                                            dom = 'Bfrtip',
                                            buttons =
                                             list('colvis', list(
                                               extend = 'collection',
                                               buttons = list(list(extend='csv',
                                                                   filename = 'hitCount'),
                                                              list(extend='excel',
                                                                   filename = 'hitCount'),
                                                              list(extend='pdf',
                                                                   filename= 'hitCount')),
                                               text = 'Download')
                                             )
                                           ))
  
  tableGroup <- formatStyle(tableGroup, names(statsOfGroupOrdered)[maxChem],
                            background = styleColorBar(range(statsOfGroupOrdered[,maxChem],na.rm=TRUE), 'goldenrod'),
                            backgroundSize = '100% 90%',
                            backgroundRepeat = 'no-repeat',
                            backgroundPosition = 'center' )
  
  tableGroup <- formatStyle(tableGroup, names(statsOfGroupOrdered)[meanChem],
                            background = styleColorBar(range(statsOfGroupOrdered[,meanChem],na.rm=TRUE), 'wheat'),
                            backgroundSize = '100% 90%',
                            backgroundRepeat = 'no-repeat',
                            backgroundPosition = 'center')

  return(tableGroup)
}


#' statsOfGroup
#' 
#' Summarize data for most graphs/tables
#' @param chemicalSummary data frame
#' @param category character
#' @param hit_threshold numeric
#' @export
#' @keywords internal
statsOfGroup <- function(chemicalSummary, category, hit_threshold){
  
  Class <- Bio_category <- site <- EAR <- sumEAR <- chnm <- n <- hits <- ".dplyr"
  
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
  
  statsOfGroup <- chemicalSummary %>%
    group_by(site, date,category) %>%
    summarise(sumEAR = sum(EAR),
              hits = sum(EAR > hit_threshold)) %>%
    group_by(site,category) %>%
    summarise(`Hits per Sample` = sum(sumEAR > hit_threshold),
              # `Mean(n) of hits` = sum(mean(sumEAR) > hit_threshold),
              `Individual Hits` = sum(hits),
              nSamples = n()) %>%
    data.frame()
  
  return(statsOfGroup)
}