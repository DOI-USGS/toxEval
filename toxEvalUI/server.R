library(DT)
library(dplyr)
library(toxEval)

shinyServer(function(input, output) {
  
#   downloaded <- function(i){
#     
#     if (i <= length(input$dataset)){
#       data = xts(NULL)
#       metadata <- get_var_src_codes(var_src==input$dataset[i], out=c('units','var_descrip','var'))
#       ylab <- paste0(metadata$var_descrip, ' (',metadata$units,')')
#       
#       for (site in input$site){
#         rdfile <- file.path(tempdir(),paste0(site,'-ts_',input$dataset[i],'.RData'))
#         if (file.exists(rdfile)){
#           message("file is cached, using local:", input$dataset[i])
#           load(file = rdfile) # loads tsobject
#         } else {
#           message("downloading file:", input$dataset[i])
#           file <- download_ts(var_src = input$dataset[i], site_name = site, on_local_exists = 'skip', on_remote_missing = "return_NA")
#           if (is.na(file))
#             tsobject <- xts(NULL)
#           else {
#             unitt <- read_ts(file)
#             var <- v(unitt)
#             tsobject <- xts(var[[2]], order.by = as.POSIXct(var[['DateTime']]), tz='UTC')
#           }
#           save(tsobject, list='tsobject', file = rdfile)
#         }
#         
#         oldnames <- names(data)
#         data <- merge(data, tsobject)
#         if(length(tsobject) > 0)
#           names(data) <- c(oldnames,paste0(site,' (',metadata$units,')'))
#         
#       }
#       list(data=data, ylab=ylab, units = metadata$units, var=metadata$var)
#     } else {
#       NULL
#     }# // else do nothing
#   }
#   
#   buildDy <- function(i){
#     ts <- downloaded(i)
#     if (!is.null(ts) && length(ts$data) > 0 ){
#       dygraph(ts$data, group = "powstreams") %>%
#         dyOptions(colors = RColorBrewer::brewer.pal(max(3,length(names(ts$data))), "Dark2")) %>%
#         dyHighlight(highlightSeriesBackgroundAlpha = 0.65, hideOnMouseOut = TRUE) %>%
#         dyAxis("y", label = ts$ylab)
#     }
#   }
#   
#   output[['dygraph1']] <- renderDygraph({
#     buildDy(1)
#   })
#   output[['dygraph2']] <- renderDygraph({
#     buildDy(2)
#   })
#   output[['dygraph3']] <- renderDygraph({
#     buildDy(3)
#   })
  output$groupControl <- renderUI({
    
    selectInput("group", label = "Group in column",
                choices = unique(endPointInfo[,input$groupCol])[!is.na(unique(endPointInfo[,input$groupCol]))],
                multiple = FALSE)
  })
})