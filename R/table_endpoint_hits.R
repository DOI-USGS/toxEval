#' table_endpoint_hits
#' 
#' Table of ranks
#' @param chemicalSummary data frame from \code{get_chemical_summary}
#' @param mean_logic logical \code{TRUE} is mean, \code{FALSE} is maximum
#' @param category either "Biological", "Chemical Class", or "Chemical"
#' @param hit_threshold numeric threshold defining a "hit"
#' @export
#' @import DT
#' @rdname table_endpoint_hits
#' @importFrom stats median
#' @importFrom tidyr spread unite
#' @importFrom dplyr full_join filter mutate select left_join right_join
#' @examples
#' # This is the example workflow:
#' path_to_tox <-  system.file("extdata", package="toxEval")
#' file_name <- "OWC_data_fromSup.xlsx"
#'
#' full_path <- file.path(path_to_tox, file_name)
#' \dontrun{
#' tox_list <- create_toxEval(full_path)
#' 
#' ACClong <- get_ACC(tox_list$chem_info$CAS)
#' ACClong <- remove_flags(ACClong)
#' 
#' cleaned_ep <- clean_endPoint_info(endPointInfo)
#' filtered_ep <- filter_groups(cleaned_ep)
#' chemicalSummary <- get_chemical_summary(tox_list, ACClong, filtered_ep)
#' 
#' hits_df <- endpoint_hits(chemicalSummary, category = "Biological")                        
#' table_endpoint_hits(chemicalSummary, category = "Biological")
#' table_endpoint_hits(chemicalSummary, category = "Chemical Class")
#' table_endpoint_hits(chemicalSummary, category = "Chemical")
#' }
table_endpoint_hits <- function(chemicalSummary, 
                           category = "Biological",
                           mean_logic = FALSE,
                           hit_threshold = 0.1){
  
  chnm <- CAS <- ".dplyr"
  
  fullData <- endpoint_hits(chemicalSummary = chemicalSummary,
                           category = category,
                           mean_logic = mean_logic,
                           hit_threshold = hit_threshold)

  if(category == "Chemical"){
    casKey <- select(chemicalSummary, chnm, CAS) %>%
      distinct()

    hits <- sapply(fullData, function(x) as.character(x))

    for(k in 1:nrow(fullData)){
      for(z in 2:ncol(fullData)){
        if(!is.na(fullData[k,z])){
          if(fullData[k,z] < 10){
            hit_char <- paste0("0",fullData[k,z])
          } else{
            hit_char <- as.character(fullData[k,z])
          }
          hits[k,z] <- paste(hit_char,createLink(cas = casKey$CAS[casKey$chnm == names(fullData)[z]],
                                  endpoint = fullData[k,1]))
        }
      }
    }

    fullData <- hits
  }
  
  fullData <- DT::datatable(fullData, extensions = 'Buttons',
                              escape = FALSE,
                              rownames = FALSE,
                              options = list(dom = 'Bfrtip',
                                             buttons = list('colvis'),
                                             scrollX = TRUE,
                                             order=list(list(1,'desc'))))
  return(fullData)
}

#' @rdname table_endpoint_hits
#' @export
endpoint_hits <- function(chemicalSummary, 
                         category = "Biological",
                         mean_logic = FALSE,
                         hit_threshold = 0.1){
  Bio_category <- Class <- EAR <- sumEAR <- value <- calc <- chnm <- choice_calc <- n <- nHits <- site <- ".dplyr"
  endPoint <- meanEAR <- nSites <- CAS <- ".dplyr"
  
  match.arg(category, c("Biological","Chemical Class","Chemical"))
  
  fullData_init <- data.frame(endPoint="",stringsAsFactors = FALSE)
  fullData <- fullData_init
  
  if(category == "Chemical"){
    chemicalSummary <- mutate(chemicalSummary, category = chnm)
  } else if (category == "Chemical Class"){
    chemicalSummary <- mutate(chemicalSummary, category = Class)
  } else {
    chemicalSummary <- mutate(chemicalSummary, category = Bio_category)
  }
  
  if(length(unique(chemicalSummary$site)) > 1){
    
    fullData <- chemicalSummary %>%
      group_by(site, category, endPoint, date) %>%
      summarize(sumEAR = sum(EAR)) %>%
      group_by(site, category, endPoint) %>%
      summarize(meanEAR = ifelse(mean_logic, mean(sumEAR),max(sumEAR))) %>%
      group_by(category, endPoint) %>%
      summarize(nSites = sum(meanEAR > hit_threshold)) %>%
      spread(category, nSites)
    
  } else {
    
    fullData <- chemicalSummary %>%
      group_by(category, endPoint, date) %>%
      summarize(sumEAR = sum(EAR)) %>%
      group_by(category, endPoint) %>%
      summarize(meanEAR = ifelse(mean_logic, mean(sumEAR),max(sumEAR))) %>%
      group_by(category, endPoint) %>%
      summarize(nSites = sum(meanEAR > hit_threshold)) %>%
      spread(category, nSites)
    
  }

  fullData <- fullData[(rowSums(fullData[,-1],na.rm = TRUE) != 0),]
  fullData <- fullData[, colSums(is.na(fullData)) != nrow(fullData)]

  sumOfColumns <- colSums(fullData[c(-1)],na.rm = TRUE)
  if(!all(sumOfColumns == 0)){
    orderData <- order(sumOfColumns,decreasing = TRUE)
    orderData <- orderData[sumOfColumns[orderData] != 0] + 1
    
    fullData <- fullData[,c(1,orderData)]   
  }
  
  fullData <- fullData[order(fullData[[2]], decreasing = TRUE),]
  
  return(fullData)
}

#' createLink
#' 
#' Create links
#' @param cas character
#' @param endpoint character
#' @param hits character
#' @export
#' @keywords internal
createLink <- function(cas, endpoint) {
  paste0('<a href="http://actor.epa.gov/dashboard/#selected/',cas,"+",endpoint,'" target="_blank">&#9432;</a>')
}