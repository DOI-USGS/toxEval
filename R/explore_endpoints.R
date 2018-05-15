#' Explore data in the Shiny Application
#' 
#' Open an interactive app in a browser. See the "Shiny App" vignette:
#' \url{../doc/shinyApp.html} for more details. Using this 
#' function is a quick and convenient way
#' to explore data. For more customization, the R-code to 
#' produce each graph and table is displayed in the app. That is 
#' a good starting-point for a custom analysis.
#' 
#' @param browse Logical. Use browser for running Shiny app.
#' @export
#' @importFrom shiny runApp
#' @importFrom tools file_ext
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
#' @import shinydashboard
#' @import shinyAce
#' @examples 
#' \dontrun{
#' explore_endpoints()
#' }
explore_endpoints <- function(browse=TRUE){
  runApp(system.file('shiny', package='toxEval'), launch.browser = browse)
}