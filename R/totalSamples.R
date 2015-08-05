#' Create total samples info
#'
#' Create total samples
#'
#' @param chemicalSummary data frame returned from \code{chemSummBasic}
#' @param chemicalSummary.Key character column name for grouping by (default = "endPoint")
#' @param chemicalSummary.site character column name for group by site
#' @export
#' @examples
#' wData <- wData
#' pCodeInfo <- pCodeInfo
#' endPoint <- endPointToxCreate(pCodeInfo)
#' chemicalSummary <- chemSummBasic(wData,pCodeInfo,endPoint)
#' ts <- totalSamples(chemicalSummary)
totalSamples <- function(chemicalSummary, chemicalSummary.Key="endPoint",chemicalSummary.site="site"){
  
  totalSamples_ep <- select_(chemicalSummary,chemicalSummary.site,chemicalSummary.Key) %>%
    distinct()%>%
    group_by_(chemicalSummary.Key) %>%
    summarise(totalSites=n())%>%
    rename_("endPoint"=chemicalSummary.Key) %>%
    mutate(endPoint=as.character(endPoint))
  
  return(totalSamples_ep)
  
}