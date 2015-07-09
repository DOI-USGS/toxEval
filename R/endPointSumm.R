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
#' @import dplyr
#' @export
#' @examples
#' wData <- wData
#' pCodeInfo <- pCodeInfo
#' endPoint <- endPointToxCreate(pCodeInfo)
#' chemicalSummary <- chemSummBasic(wData,pCodeInfo,endPoint)
#' endPointInfo <- endPointInfo
#' epSumm <- endPointSumm(chemicalSummary, endPointInfo=endPointInfo)
endPointSumm <- function(chemicalSummary, chemicalSummary.Key="endPoint",chemicalSummary.site="site",
                         endPointInfo, endPointInfo.Key = "assay_component_endpoint_name",
                         EAR.key="EAR"){
  
  totalSamples_ep <- select_(chemicalSummary,chemicalSummary.site,chemicalSummary.Key) %>%
    distinct()%>%
    group_by_(chemicalSummary.Key) %>%
    summarise(totalSites=n())%>%
    rename_("endPoint"=chemicalSummary.Key) %>%
    mutate(endPoint=as.character(endPoint))
  
  
  endpointSummary <- chemicalSummary %>%
    mutate_("hits"= paste0("as.numeric(",EAR.key," > 0.1)")) %>%
    rename_("EAR"=EAR.key) %>%
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

