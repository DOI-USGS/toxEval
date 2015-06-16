#' Create site summary
#'
#' Create site summary
#'
#' @param chemicalSummary data frame returned from \code{chemSummBasic}
#' @param newSiteKey named vector 
#' @import  dplyr
#' @export
#' @examples
#' wData <- wData
#' pCodeInfo <- pCodeInfo
#' packagePath <- system.file("extdata", package="toxEval")
#' filePath <- file.path(packagePath, "stationINFO.RData")
#' load(file=filePath)
#' newSiteKey <- setNames(stationINFO$shortName, stationINFO$fullSiteID)
#' endPoint <- endPointToxCreate(pCodeInfo)
#' chemicalSummary <- chemSummBasic(wData,pCodeInfo,endPoint)
#' siteSummary <- siteSumm(chemicalSummary,newSiteKey)
siteSumm <- function(chemicalSummary,newSiteKey){

  siteSummary <- chemicalSummary %>%
    mutate(site=as.character(newSiteKey[site])) %>% 
    select(site, chnm, EAR, endPoint, class, date) %>% 
    arrange(site, chnm, EAR)
  
  summ1 <- siteSummary %>%
    group_by(site, date) %>%  #This ignore chemicals...needs to be here for freq
    summarize(hits= as.numeric(any(EAR > 0.1))) %>%
    group_by(site) %>%
    summarize(freq=sum(hits)/n_distinct(date),
              nSamples=n_distinct(date)) 
  
  summ2 <-  siteSummary %>%
    group_by(site) %>%
    summarize(nChem = length(unique(chnm[EAR > 0.1])),
              nEndPoints = length(unique(endPoint[EAR > 0.1])),
              maxEAR = max(EAR))
  
  summary <- left_join(summ1, summ2, by="site") %>%
    arrange(desc(maxEAR))
  
  return(summary)
}