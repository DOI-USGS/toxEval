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
#' @importFrom dplyr n
#' @importFrom grDevices colorRampPalette
#' @importFrom leaflet colorBin
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
#' mapData <- getMapInfo(chemicalSummary, tox_list$chem_site, "Biological") 
#' }
getMapInfo <- function(chemicalSummary,
                    chem_site,
                    category = "Biological",
                    mean_logic = FALSE){

  match.arg(category, c("Biological","Chemical Class","Chemical"))
  
  site <- meanEAR <- nSamples <- `Short Name` <- dec_lat <- dec_lon <- ".dplyr"
  
  siteToFind <- unique(chemicalSummary$shortName)
  
  if(category == "Biological"){
    typeWords <- "groups"
  } else if (category == "Chemical"){
    typeWords <- "chemicals"
  } else {
    typeWords <- "chemical classes"
  }
  
  mapData <- chem_site[chem_site$`Short Name` %in% siteToFind,
                       c("Short Name", "dec_lat", "dec_lon", "SiteID")]
  
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

#' makeMap
#' 
#' makeMap
#' 
#' @param chemicalSummary data frame from \code{get_chemical_summary}
#' @param category either "Biological", "Chemical Class", or "Chemical"
#' @param mean_logic logical \code{TRUE} is mean, \code{FALSE} is maximum
#' @param chem_site data frame with at least columns SiteID, site_grouping, and Short Name
#' @export
#' @import leaflet
#' @examples
#' # This is the example workflow:
#' path_to_tox <-  system.file("extdata", package="toxEval")
#' file_name <- "OWC_data_fromSup.xlsx"
#'
#' full_path <- file.path(path_to_tox, file_name)
#' tox_list <- create_toxEval(full_path)
#' 
#' \dontrun{
#' 
#' ACClong <- get_ACC(tox_list$chem_info$CAS)
#' ACClong <- remove_flags(ACClong)
#' 
#' cleaned_ep <- clean_endPoint_info(endPointInfo)
#' filtered_ep <- filter_groups(cleaned_ep)
#' chemicalSummary <- get_chemical_summary(tox_list, ACClong, filtered_ep)
#' 
#' makeMap(chemicalSummary, tox_list$chem_site, "Biological")   
#' makeMap(chemicalSummary, tox_list$chem_site, "Chemical Class")
#' makeMap(chemicalSummary, tox_list$chem_site, "Chemical") 
#' }
makeMap <- function(chemicalSummary,
                    chem_site,
                    category = "Biological",
                    mean_logic = FALSE){
  
  maxEARWords <- ifelse(mean_logic,"meanEAR","maxEAR")
  
  mapDataList <- getMapInfo(chemicalSummary, 
                            chem_site = chem_site, 
                            category = category,
                            mean_logic = mean_logic)
  
  mapData <- mapDataList$mapData
  pal <- mapDataList$pal
  siteToFind <- unique(chemicalSummary$site)
  
  if(length(siteToFind) == 1){
    
    mapData <- filter(chem_site, SiteID == siteToFind) %>%
      mutate(nSamples = median(mapData$count),
             meanMax = median(mapData$meanMax),
             sizes = median(mapData$sizes))
  }
  map <- leaflet::leaflet(height = "500px", data=mapData) %>%
    leaflet::addProviderTiles("CartoDB.Positron") %>%
    leaflet::setView(lng = -83.5, lat = 44.5, zoom=6) %>%
    leaflet::clearMarkers() %>%
    leaflet::clearControls() %>%
    leaflet::setView(lng = mean(mapData$dec_lon, na.rm = TRUE), 
            lat = mean(mapData$dec_lat, na.rm = TRUE), zoom=6) %>%
    leaflet::fitBounds(lng1 = min(mapData$dec_lon, na.rm = TRUE), 
              lat1 = min(mapData$dec_lat, na.rm = TRUE), 
              lng2 = max(mapData$dec_lon, na.rm = TRUE), 
              lat2 = max(mapData$dec_lat, na.rm = TRUE)) %>%
    leaflet::addCircleMarkers(lat=~dec_lat, lng=~dec_lon,
                     popup=paste0('<b>',mapData$`Short Name`,"</b><br/><table>",
                                  "<tr><td>",maxEARWords,": </td><td>",sprintf("%.1f",mapData$meanMax),'</td></tr>',
                                  "<tr><td>Number of Samples: </td><td>",mapData$count,'</td></tr>',
                                  '</table>') ,
                     fillColor = ~pal(meanMax),
                     fillOpacity = 0.8,
                     radius = ~sizes,
                     stroke=FALSE,
                     opacity = 0.8)
  
  if(length(siteToFind) > 1){
    map <- leaflet::addLegend(map,pal = pal,
                     position = 'bottomleft',
                     values=~meanMax,
                     opacity = 0.8,
                     labFormat = leaflet::labelFormat(digits = 2), #transform = function(x) as.integer(x)),
                     title = paste("Sum of",category,"EAR<br>",
                                   ifelse(mean_logic,'Mean','Max'),"at site"))
  }
  
  return(map)
  
}
