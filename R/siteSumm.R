#' Create site summary
#'
#' Create site summary
#'
#' @param wData data frame with a date column (can be NA) 
#' whose name is defined by the date argument, station column whose name is
#' defined by the station argument, and the remaining columns are measurements that
#' contain a chemKey (such as parameter code)
#' @param pCodeInfoDF data frame
#' @param endPoint data frame 
#' @param date character name of date column in wData
#' @param station character name of station column in wData
#' @param casrn_pCode character name of casrn column in pCodeInfo
#' @param class_pCode character name of class column in pCodeInfo
#' @param code_pCode character name of code column in pCodeInfo
#' @param casrn_ep character name of class column in endPoint
#' @import dplyr
#' @import tidyr
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