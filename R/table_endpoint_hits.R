#' table_endpoint_hits
#' 
#' Table of ranks
#' @param chemicalSummary data frame from \code{get_chemical_summary}
#' @param mean_logic logical \code{TRUE} is mean, \code{FALSE} is maximum
#' @param category either "Biological", "Chemical Class", or "Chemical"
#' @param hit_threshold numeric threshold defining a "hit"
#' @export
#' @import DT
#' @importFrom stats median
#' @importFrom tidyr spread unite
#' @importFrom dplyr full_join filter mutate select left_join right_join
#' @examples
#' library(readxl)
#' path_to_tox <-  system.file("extdata", package="toxEval")
#' file_name <- "OWC_data_fromSup.xlsx"
#' full_path <- file.path(path_to_tox, file_name)
#' 
#' chem_data <- read_excel(full_path, sheet = "Data")
#' chem_info <- read_excel(full_path, sheet = "Chemicals") 
#' chem_site <- read_excel(full_path, sheet = "Sites")
#' ACClong <- get_ACC(chem_info$CAS)
#' ACClong <- remove_flags(ACClong)
#' 
#' cleaned_ep <- clean_endPoint_info(endPointInfo)
#' filtered_ep <- filter_groups(cleaned_ep)
#' 
#' chemicalSummary <- get_chemical_summary(ACClong,
#'                                         filtered_ep,
#'                                        chem_data, 
#'                                         chem_site, 
#'                                         chem_info)
#' table_endpoint_hits(chemicalSummary, category = "Biological")
#' table_endpoint_hits(chemicalSummary, category = "Chemical Class")
#' table_endpoint_hits(chemicalSummary, category = "Chemical")
table_endpoint_hits <- function(chemicalSummary, 
                           category = "Biological",
                           mean_logic = FALSE,
                           hit_threshold = 0.1){
  
  Bio_category <- Class <- EAR <- sumEAR <- value <- calc <- chnm <- choice_calc <- n <- nHits <- site <- ".dplyr"
  endPoint <- meanEAR <- nSites <- casrn <- ".dplyr"
  
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
  
    for(i in unique(chemicalSummary$category)){
      
      dataSub <- chemicalSummary %>%
        filter(category == i) %>%
        group_by(site, category, endPoint, date) %>%
        summarize(sumEAR = sum(EAR)) %>%
        group_by(site, category, endPoint) %>%
        summarize(meanEAR = ifelse(mean_logic, mean(sumEAR),max(sumEAR))) %>%
        group_by(category, endPoint) %>%
        summarize(nSites = sum(meanEAR > hit_threshold)) %>%
        data.frame() %>%
        arrange(desc(nSites)) %>%
        spread(category, nSites)

      if(ncol(dataSub) > 2){
        dataSub <- dataSub[,c(1,1+which(colSums(dataSub[,-1],na.rm = TRUE) != 0))]
      }
  
      if(is.data.frame(dataSub)){
        if(ncol(dataSub) > 2){
          dataSub <- dataSub[(rowSums(dataSub[,-1],na.rm = TRUE) != 0),]
        } else {
          dataSub <- dataSub[which(dataSub[,-1] != 0 ),]
        }
        fullData <- full_join(fullData,dataSub, by="endPoint")
      }
    }
  
  } else {
  
    for(i in unique(chemicalSummary$category)){
      dataSub <- chemicalSummary %>%
        filter(category == i) %>%
        group_by(category, endPoint, date) %>%
        summarise(sumEAR=sum(EAR)) %>%
        data.frame() %>%
        group_by(endPoint, category) %>%
        summarise(nSites = sum(sumEAR > hit_threshold))%>%
        data.frame() %>%
        arrange(desc(nSites)) %>%
        spread(category, nSites)

      if(ncol(dataSub) > 2){
        dataSub <- dataSub[,c(1,1+which(colSums(dataSub[,-1],na.rm = TRUE) != 0))]
      }
  
      if(is.data.frame(dataSub)){
        if(ncol(dataSub) > 2){
          dataSub <- dataSub[(rowSums(dataSub[,-1],na.rm = TRUE) != 0),]
        } else {
          dataSub <- dataSub[which(dataSub[,-1] != 0 ),]
        }

        fullData <- full_join(fullData,dataSub, by="endPoint")
      }
    }
  }

  sumOfColumns <- colSums(fullData[c(-1)],na.rm = TRUE)
  if(!all(sumOfColumns == 0)){
    orderData <- order(sumOfColumns,decreasing = TRUE)
    orderData <- orderData[sumOfColumns[orderData] != 0] + 1
    
    fullData <- fullData[,c(1,orderData)]   
  }

  if(category == "Chemical"){
    casKey <- select(chemicalSummary, chnm, casrn) %>%
      distinct()

    hits <- sapply(fullData, function(x) as.character(x))

    for(k in 1:nrow(fullData)){
      for(z in 2:ncol(fullData)){
        if(!is.na(fullData[k,z])){
          hits[k,z] <- createLink(cas = casKey$casrn[casKey$chnm == names(fullData)[z]],
                                  endpoint = fullData[k,1],
                                  hits = fullData[k,z])
        }
      }
    }

    fullData <- hits
  }
  
  fullData <- DT::datatable(fullData, extensions = 'Buttons',
                              escape = FALSE,
                              rownames = FALSE,
                              options = list(dom = 'Bfrtip',
                                             buttons =
                                               list('colvis', list(
                                                 extend = 'collection',
                                                 buttons = list(list(extend='csv',
                                                                     filename = 'epHits'),
                                                                list(extend='excel',
                                                                     filename = 'epHits'),
                                                                list(extend='pdf',
                                                                     filename= 'epHits')),
                                                 text = 'Download'
                                               )),
                                             scrollX = TRUE,
                                             pageLength = nrow(fullData),
                                             order=list(list(2,'desc'))))
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
createLink <- function(cas, endpoint, hits) {
  paste0('<a href="http://actor.epa.gov/dashboard/#selected/',cas,"+",endpoint,'" target="_blank" >',hits,'</a>')
}