#' Explore endpoint groupings
#' 
#' Open an interactive app
#' 
#' @param browse use browser for map rendering
#' @export
#' @importFrom shiny runApp
#' @importFrom leaflet leaflet
#' @importFrom leaflet addLegend
#' @importFrom leaflet addCircles
#' @importFrom leaflet clearControls
#' @importFrom leaflet clearShapes
#' @import DT
#' @import ggplot2
#' @import data.table
#' @import RColorBrewer
#' @import grid
#' @import gridExtra 
explore_endpoints <- function(browse=TRUE){
  runApp(system.file('shiny', package='toxEval'), launch.browser = browse)
}