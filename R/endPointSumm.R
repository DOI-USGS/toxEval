#' Create chemical summary
#'
#' Create chemical summary
#'
#' @param chemicalSummary data frame returned from \code{chemSummBasic}
#' @param chemicalSummary.Key character column name for grouping by (default = "endPoint")
#' @param chemicalSummary.site character column name for group by site
#' @param endPointInfo dataframe
#' @param endPointInfo.Key column name of endpoint
#' @param EAR.key column name in chemicalSummary of EAR
#' @export
#' @examples
#' \dontrun{
#' wData <- wData
#' pCodeInfo <- pCodeInfo
#' endPoint <- endPointToxCreate(pCodeInfo)
#' chemicalSummary <- chemSummBasic(wData,pCodeInfo,endPoint)
#' endPointInfo <- endPointInfo
#' epSumm <- endPointSumm(chemicalSummary, endPointInfo=endPointInfo)
#' }
endPointSumm <- function(chemicalSummary, chemicalSummary.Key="endPoint",chemicalSummary.site="site",
                         endPointInfo, endPointInfo.Key = "assay_component_endpoint_name",
                         EAR.key="EAR"){
  
  totalSamples_ep <- totalSamples(chemicalSummary = chemicalSummary, chemicalSummary.Key=chemicalSummary.Key,
                                  chemicalSummary.site=chemicalSummary.site)
  
  
  endpointSummary <- chemicalSummary %>%
    group_by_(chemicalSummary.site, chemicalSummary.Key) %>%
    summarize(hits=as.numeric(any(hits > 0)),
              maxEAR=max(EAR)) %>%
    group_by_(chemicalSummary.Key) %>%
    summarize(nSites=sum(hits),
              maxEAR=max(maxEAR)) %>%
    rename_("endPoint"=chemicalSummary.Key) %>%
    mutate(endPoint=as.character(endPoint))%>%
    left_join(endPointInfo, by=c("endPoint"=endPointInfo.Key)) %>% 
    left_join(totalSamples_ep, by=c("endPoint")) %>%
    mutate(freq=nSites/totalSites)%>%
    data.frame %>%
    arrange(desc(freq)) 
  
  return(endpointSummary)
}

