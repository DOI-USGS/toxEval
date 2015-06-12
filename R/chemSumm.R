#' Create chemical summary
#'
#' Create chemical summary
#'
#' @param wData data frame with a date column (can be NA) 
#' whose name is defined by the date argument, station column whose name is
#' defined by the station argument, and the remaining columns are measurements that
#' contain a chemKey (such as parameter code)
#' @param pCodeInfo data frame
#' @param endPoint data frame 
#' @param date character name of date column in wData
#' @param station character name of station column in wData
#' @param casrn_pCode character name of casrn column in pCodeInfo
#' @param class_pCode character name of class column in pCodeInfo
#' @param code_pCode character name of code column in pCodeInfo
#' @param casrn_ep character name of class column in endPoint
#' @return chemicalSummary data frame
#' @import dplyr
#' @import tidyr
#' @export
#' @examples
#' wData <- wData
#' pCodeInfo <- pCodeInfo
#' 
chemSumm <- function(wData, pCodeInfo,endPoint,
                     date="ActivityStartDateGiven", station="site",
                     casrn_pCode="casrn",class_pCode="class",code_pCode="parameter_cd",
                     casrn_ep="casn"){
  
  chemicalSummary <- wData %>%
    gather_("pCode", "measuredValue", gather_cols=names(wData)[!(names(wData) %in% c(date,station))])  %>%
    rename_("date"=date) %>%
    filter(!is.na(measuredValue)) %>%
    separate(pCode, into=c("colHeader","pCode"), "_") %>% #totally not generalized
    left_join(pCodeInfo[c(code_pCode,casrn_pCode,class_pCode)], by=c("pCode"=code_pCode))  %>%
    right_join(endPoint, by=c(casrn=casrn_ep)) %>%
    select(-mlWt, -conversion, -casrn,  -Units, -pCode, -colHeader) %>%
    gather(endPoint, endPointValue, -class, -site, -measuredValue, -chnm, -date)  %>% 
    filter(!is.na(endPointValue)) %>%
    mutate(EAR=measuredValue/endPointValue)
  
  
    chemSum1 <- chemicalSummary %>%
      mutate(hits= as.numeric(EAR > 0.1)) %>%
      group_by(chnm,  class, site) %>%
      summarize(hits=as.numeric(any(hits > 0)), 
                median=median(EAR),
                maxEAR=max(EAR)) %>%
      group_by(chnm, class) %>%
      summarize(freq=sum(hits)/n_distinct(site),
                nSites=sum(hits),
                aveMedian=mean(median)) 
    
    totalSamples <- select(chemicalSummary,date,site,chnm) %>%
      distinct()%>%
      group_by(chnm) %>%
      summarise(nSamples=n()) %>%
      mutate(chnm=as.character(chnm))
    
    chemSum2 <- chemicalSummary %>%
      group_by(chnm,  class) %>%  
      summarize(maxEAR=max(EAR),
                nEndPoints=length(unique(endPoint[EAR > 0.1])),
                EAR90overall=quantile(EAR,probs = .9)) %>%
      left_join(chemSum1, by=c("chnm","class")) %>%
      left_join(totalSamples, by=c("chnm")) %>%
      data.frame %>%
      arrange(desc(freq)) %>%
      select(chnm, class, freq, maxEAR, nEndPoints, nSites, nSamples)
  
  
  return(chemSum2)
  
}


  