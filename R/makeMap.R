#' getMapInfo
#' 
#' getMapInfo
#' 
#' @param chemicalSummary data frame from \code{get_chemical_summary}
#' @param category either "Biological", "Chemical Class", or "Chemical"
#' @param chem_site data frame with at least columns SiteID, site_grouping, and Short Name
#' @param mean_logic logical \code{TRUE} is mean, \code{FALSE} is maximum
#' @export
#' @import ggplot2
#' @importFrom stats median
#' @importFrom grDevices colorRampPalette
#' @importFrom leaflet colorBin
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
#' mapData <- getMapInfo(chemicalSummary, chem_site, "Biological") 
getMapInfo <- function(chemicalSummary,
                    chem_site,
                    category = "Biological",
                    mean_logic = FALSE){

  match.arg(category, c("Biological","Chemical Class","Chemical"))
  
  site <- nSamples <- `Short Name` <- Fullname <- dec_lat <- dec_lon <- ".dplyr"
  
  siteToFind <- chem_site$`Short Name`
  
  if(category == "Biological"){
    typeWords <- "groups"
  } else if (category == "Chemical"){
    typeWords <- "chemicals"
  } else {
    typeWords <- "chemical classes"
  }
  
  mapData <- chem_site[,c("Short Name", "dec_lat", "dec_lon", "SiteID")]
  
  nSamples <- select(chemicalSummary,site,date) %>%
    distinct() %>%
    group_by(site) %>%
    summarize(count = n())
  
  meanStuff <- graphData(chemicalSummary = chemicalSummary, 
            category = category, mean_logic = mean_logic) %>%
    group_by(site) %>%
    summarize(meanMax = max(meanEAR)) %>%
    left_join(nSamples, by="site")
  
  mapData <- left_join(mapData, meanStuff, by=c("SiteID"="site"))
  
  col_types <- c("darkblue","dodgerblue","green4","gold1","orange","brown","red")

  counts <- mapData$count       
  
  if(length(siteToFind) > 1){
    leg_vals <- unique(as.numeric(quantile(mapData$meanMax, probs=c(0,0.01,0.1,0.25,0.5,0.75,0.9,.99,1), na.rm=TRUE)))
    pal = colorBin(col_types, mapData$meanMax, bins = leg_vals)
    rad <-3*seq(1,4,length.out = 16)
    
    if(sum(mapData$count, na.rm = TRUE) == 0){
      mapData$sizes <- rad[1]
    } else {
      mapData$sizes <- rad[as.numeric(cut(mapData$count, breaks=16))]
    }
    
  } else {
    leg_vals <- unique(as.numeric(quantile(c(0,mapData$meanMax), probs=c(0,0.01,0.1,0.25,0.5,0.75,0.9,.99,1), na.rm=TRUE)))
    pal = colorBin(col_types, c(0,mapData$meanMax), bins = leg_vals)
    mapData$sizes <- 3
  }
  
  return(list(mapData=mapData, pal=pal))
  
}
