#' Create chemical summary
#'
#' Create chemical summary
#'
#' @param chemicalSummary data frame returned from \code{chemSummBasic}
#' @param EAR.key column name of EAR
#' @param chnmCol column name of chemical names
#' @param classCol column name of chemical class
#' @param siteCol column name of site
#' @param dateCol column name of date
#' @export
#' @examples
#' wData <- wData
#' pCodeInfo <- pCodeInfo
#' endPoint <- endPointToxCreate(pCodeInfo)
#' chemicalSummary <- chemSummBasic(wData,pCodeInfo,endPoint)
#' chemSum1 <- chemSumm(chemicalSummary)
chemSumm <- function(chemicalSummary,EAR.key="EAR",chnmCol="chnm",
                     classCol="class",siteCol="site",dateCol="date"){

  chemSum1 <- chemicalSummary %>%
    mutate_("hits"= paste0("as.numeric(",EAR.key," > 0.1)")) %>%
    rename_("EAR"=EAR.key, "site"=siteCol,"class"=classCol,"chnm"=chnmCol) %>%
    group_by(chnm,  class, site) %>%
    summarize(hits=as.numeric(any(hits > 0)),
              median=median(EAR),
              maxEAR=max(EAR)) %>%
    group_by(chnm, class) %>%
    summarize(freq=sum(hits)/n_distinct(site),
              nSites=sum(hits),
              aveMedian=mean(median)) 
  
  totalSamples <- select_(chemicalSummary,dateCol,siteCol,chnmCol) %>%
    rename_("chnm"=chnmCol) %>%
    distinct()%>%
    group_by(chnm) %>%
    summarise(nSamples=n()) %>%
    mutate(chnm=as.character(chnm))
  
  chemSum2 <- chemicalSummary %>%
    rename_("chnm"=chnmCol, "class"=classCol, "EAR"=EAR.key) %>%
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


#' Create basic chemical summary
#'
#' Create basic chemical summary
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
#' @return chemicalSummary data frame
#' @importFrom tidyr gather_ 
#' @importFrom tidyr gather
#' @importFrom tidyr separate
#' @export
#' @examples
#' wData <- wData
#' pCodeInfo <- pCodeInfo
#' endPoint <- endPointToxCreate(pCodeInfo)
#' chemicalSummary <- chemSummBasic(wData,pCodeInfo,endPoint)
chemSummBasic <- function(wData, pCodeInfoDF,endPoint,
                     date="ActivityStartDateGiven", station="site",
                     casrn_pCode="casrn",class_pCode="class",code_pCode="parameter_cd",
                     casrn_ep="casn"){
  
  chemicalSummary <- wData %>%
    gather_("pCode", "measuredValue", gather_cols=names(wData)[!(names(wData) %in% c(date,station))]) %>%
    rename_("date"=date) %>%
    filter(!is.na(measuredValue)) %>%
    separate(pCode, into=c("colHeader","pCode"), "_") %>% #totally not generalized
    left_join(pCodeInfoDF[c(code_pCode,casrn_pCode,class_pCode)], by=c("pCode"=code_pCode))  %>%
    right_join(endPoint, by=c(casrn=casrn_ep)) %>%
    select(-mlWt, -conversion, -Units, -pCode, -colHeader) %>%
    gather(endPoint, endPointValue, -class, -site, -measuredValue, -chnm, -date, -casrn)  %>% 
    filter(!is.na(endPointValue)) %>%
    mutate(EAR=measuredValue/endPointValue) %>%
    mutate(hits = as.numeric(EAR > 0.1),
           endPoint=as.character(endPoint)) %>%
    filter(!is.na(EAR))
  return(chemicalSummary)
}

#' Create endpoint dataframe
#' 
#' Create endpoint dataframe
#'
#' @param pCodeInfoDF data frame
#' @param casrn_pCode character name of casrn column in pCodeInfo
#' @param units_pCode character name of units column in pCodeInfo
#' @param mlWt_pCode character name of molecular weight column in pCodeInfo
#' @return endPoint data frame
#' @export
#' @examples
#' pCodeInfo <- pCodeInfo
#' endPoint <- endPointToxCreate(pCodeInfo)
endPointToxCreate <- function(pCodeInfoDF,
                              casrn_pCode="casrn", units_pCode="parameter_units", 
                              mlWt_pCode="mlWt"){
  
  unitConversion <- setNames(c(10^6, 10^3, 1), c("pg/l", "ng/l", "ug/l") )
  pCodeInfoDF[,mlWt_pCode] <- pCodeInfoDF[,mlWt_pCode] * unitConversion[tolower(pCodeInfoDF[,units_pCode])]
  
  
  AC50gainTemp <- AC50gain
  AC50 <- left_join(AC50gainTemp, 
                    pCodeInfoDF[,c(casrn_pCode, units_pCode, mlWt_pCode)],
                    by= c("casn"=casrn_pCode)) %>%
    filter_(!is.na(units_pCode)) %>%
    rename_("desiredUnits" = units_pCode) %>%
    mutate_("conversion" = mlWt_pCode) %>%
    select_("casn", "chnm", "desiredUnits", mlWt_pCode, "conversion")
  
  AC50Converted <- left_join(AC50, AC50gain, by = c("casn", "chnm"))
  infoColumns <- c("casn", "chnm", "desiredUnits","mlWt", "conversion", "code","chid")
  
  endPointData <- AC50Converted[,!(names(AC50Converted) %in% infoColumns)]
  endPointData <- 10^endPointData
  endPointData <- endPointData * AC50Converted$conversion
  endPoint <- cbind(AC50, data.frame(endPointData))
  endPoint <- rename_(endPoint, "Units"="desiredUnits")
  return(endPoint)
}
  