#' Create an interactive map of the data
#' 
#' The function \code{make_tox_map} creates a \code{\link[leaflet]{leaflet}} map 
#' of the sites. This function places symbols at the location of each site in the data 
#' file that represent the magnitude of EAR (color) and the number of 
#' samples in the data set (size). This is the only function that requires 
#' "dec_lon" and "dec_lat" (decimal longitude and decimal latitude) in the 
#' data frame specified for the chem_site argument.
#' 
#' The function \code{map_tox_data} calculates the statistics for the map. It
#' my be useful on it's own.
#' 
#' @param chemical_summary Data frame from \code{\link{get_chemical_summary}}.
#' @param category Character. Either "Biological", "Chemical Class", or "Chemical".
#' @param mean_logic Logical.  \code{TRUE} displays the mean EAR from each site,
#' \code{FALSE} displays the maximum EAR from each site.
#' @param sum_logic Logical. \code{TRUE} sums the EARs in a specified grouping,
#' \code{FALSE} does not. \code{FALSE} may be better for traditional benchmarks as
#' opposed to ToxCast benchmarks.
#' @param chem_site Data frame containing the columns SiteID, site_grouping, 
#' Short Name, dec_lon, and dec_lat.
#' @export
#' @rdname make_tox_map
#' @importFrom stats quantile
#' @examples
#' # This is the example workflow:
#' path_to_tox <-  system.file("extdata", package="toxEval")
#' file_name <- "OWC_data_fromSup.xlsx"
#'
#' full_path <- file.path(path_to_tox, file_name)
#' tox_list <- create_toxEval(full_path)
#' 
#' ACC <- get_ACC(tox_list$chem_info$CAS)
#' ACC <- remove_flags(ACC)
#' 
#' cleaned_ep <- clean_endPoint_info(end_point_info)
#' filtered_ep <- filter_groups(cleaned_ep)
#' chemical_summary <- get_chemical_summary(tox_list, ACC, filtered_ep)
#' 
#' make_tox_map(chemical_summary, tox_list$chem_site, "Biological")   
#' make_tox_map(chemical_summary, tox_list$chem_site, "Chemical Class")
#' make_tox_map(chemical_summary, tox_list$chem_site, "Chemical") 
#' 
make_tox_map <- function(chemical_summary,
                    chem_site,
                    category = "Biological",
                    mean_logic = FALSE,
                    sum_logic = TRUE){
  
  SiteID <- ".dplyr"

  maxEARWords <- ifelse(mean_logic,"meanEAR","maxEAR")
  
  mapDataList <- map_tox_data(chemical_summary, 
                            chem_site = chem_site, 
                            category = category,
                            mean_logic = mean_logic,
                            sum_logic = sum_logic)
  
  mapData <- mapDataList$mapData
  pal <- mapDataList$pal
  siteToFind <- unique(chemical_summary$site)
  
  if(length(siteToFind) == 1){
    
    mapData <- dplyr::filter(chem_site, SiteID == siteToFind) %>%
      dplyr::mutate(nSamples = median(mapData$count),
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
  
  title_words <- ifelse(mean_logic,"Mean","Max")
  
  if(length(siteToFind) > 1){
    map <- leaflet::addLegend(map,pal = pal,
                     position = 'bottomleft',
                     values=~meanMax,
                     opacity = 0.8,
                     labFormat = leaflet::labelFormat(digits = 2), #transform = function(x) as.integer(x)),
                     title = paste("Sum of",category,"EAR<br>",
                                   title_words,"at site"))
  }
  
  return(map)
  
}

#' @export
#' @rdname make_tox_map
map_tox_data <- function(chemical_summary,
                       chem_site,
                       category = "Biological",
                       mean_logic = FALSE,
                       sum_logic = TRUE){
  
  match.arg(category, c("Biological","Chemical Class","Chemical"))

  site <- meanEAR <- nSamples <- `Short Name` <- dec_lat <- dec_lon <- n <- ".dplyr"
  
  siteToFind <- unique(chemical_summary$shortName)
  
  if(category == "Biological"){
    typeWords <- "groups"
  } else if (category == "Chemical"){
    typeWords <- "chemicals"
  } else {
    typeWords <- "chemical classes"
  }
  
  if(!all(c("Short Name", "dec_lat", "dec_lon", "SiteID") %in% names(chem_site))){
    stop('Map functions require columns: "Short Name", "dec_lat", "dec_lon", "SiteID" in chem_site')
  }
  
  mapData <- chem_site[chem_site$`Short Name` %in% siteToFind,
                       c("Short Name", "dec_lat", "dec_lon", "SiteID")]
  
  nSamples <- dplyr::select(chemical_summary,site,date) %>%
    dplyr::distinct() %>%
    dplyr::group_by(site) %>%
    dplyr::summarize(count = n())
  
  meanStuff <- tox_boxplot_data(chemical_summary = chemical_summary, 
                         category = category,
                         mean_logic = mean_logic,
                         sum_logic = sum_logic) %>%
    dplyr::group_by(site) %>%
    dplyr::summarize(meanMax = max(meanEAR)) %>%
    dplyr::left_join(nSamples, by="site")
  
  mapData <- dplyr::left_join(mapData, meanStuff, by=c("SiteID"="site"))
  
  col_types <- c("darkblue","dodgerblue","green4","gold1","orange","brown","red")
  
  counts <- mapData$count       
  
  if(length(siteToFind) > 1){
    leg_vals <- unique(as.numeric(quantile(mapData$meanMax, probs=c(0,0.01,0.1,0.25,0.5,0.75,0.9,.99,1), na.rm=TRUE)))
    pal = leaflet::colorBin(col_types, mapData$meanMax, bins = leg_vals)
    rad <-3*seq(1,4,length.out = 16)
    
    if(sum(mapData$count, na.rm = TRUE) == 0){
      mapData$sizes <- rad[1]
    } else {
      mapData$sizes <- rad[as.numeric(cut(mapData$count, breaks=16))]
    }
    
  } else {
    leg_vals <- unique(as.numeric(quantile(c(0,mapData$meanMax), probs=c(0,0.01,0.1,0.25,0.5,0.75,0.9,.99,1), na.rm=TRUE)))
    pal = leaflet::colorBin(col_types, c(0,mapData$meanMax), bins = leg_vals)
    mapData$sizes <- 3
  }
  
  return(list(mapData=mapData, pal=pal))
  
}