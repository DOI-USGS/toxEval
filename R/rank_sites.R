#' Rank sites by EAR
#' 
#' These functions create a table (either data frame or DT) with one row per site. 
#' There are columns of the max or mean EAR (depending on the 'mean_logic' argument) 
#' for each category ("Chemical Class", "Chemical", or "Biological"). Additionally, 
#' there are columns of the frequency of the max or mean EAR being above a user 
#' specified 'hit_threshold'.
#' 
#' The tables show slightly different results for a single site. Instead of multiple 
#' columns for category, there is now 1 row per category (since the site is known).
#' 
#' @param chemicalSummary data frame from \code{get_chemical_summary}
#' @param mean_logic character. Options are "mean", "max", or "noSum". 
#' TRUE will default to "mean" and FALSE to "max". The default value is "mean". 
#' The most appropriate use of "noSum" is for non-ToxCast benchmarks. In this case
#' the values plotted are the overall max of the sample (not the max of the sum
#' of the sample).
#' @param category either "Biological", "Chemical Class", or "Chemical"
#' @param hit_threshold numeric threshold defining a "hit"
#' @export
#' @import DT
#' @rdname rank_sites_DT
#' @importFrom stats median
#' @importFrom tidyr spread unite
#' @importFrom dplyr full_join filter mutate select left_join right_join
#' @examples
#' # This is the example workflow:
#' path_to_tox <-  system.file("extdata", package="toxEval")
#' file_name <- "OWC_data_fromSup.xlsx"
#'
#' full_path <- file.path(path_to_tox, file_name)
#' 
#' tox_list <- create_toxEval(full_path)
#' \dontrun{
#' ACClong <- get_ACC(tox_list$chem_info$CAS)
#' ACClong <- remove_flags(ACClong)
#' 
#' cleaned_ep <- clean_endPoint_info(endPointInfo)
#' filtered_ep <- filter_groups(cleaned_ep)
#' chemicalSummary <- get_chemical_summary(tox_list, ACClong, filtered_ep)
#'
#' stats_df <- rank_sites(chemicalSummary, "Biological")
#'
#' rank_sites_DT(chemicalSummary, category = "Biological")
#' rank_sites_DT(chemicalSummary, category = "Chemical Class")
#' rank_sites_DT(chemicalSummary, category = "Chemical")
#' }
rank_sites_DT <- function(chemicalSummary, 
                          category = "Biological",
                          mean_logic = FALSE,
                          hit_threshold = 0.1){
  
  Bio_category <- Class <- EAR <- maxEAR <- sumEAR <- value <- calc <- chnm <- choice_calc <- n <- nHits <- site <- ".dplyr"
  
  match.arg(category, c("Biological","Chemical Class","Chemical"))
  mean_logic <- as.character(mean_logic)
  match.arg(mean_logic, c("mean","max","noSum","TRUE","FALSE"))
  
  statsOfColumn <- rank_sites(chemicalSummary=chemicalSummary,
                                   category = category,
                                   hit_threshold = hit_threshold,
                                   mean_logic = mean_logic)

  colToSort <- 1
  
  if(mean_logic %in% c("TRUE","mean")){
    mean_logic <- TRUE
    maxEARS <- grep("meanEAR",names(statsOfColumn))
  }
  if(mean_logic %in% c("FALSE","max")){
    mean_logic <- FALSE
    maxEARS <- grep("maxEAR",names(statsOfColumn))
  }
  if(mean_logic == "noSum"){
    maxEARS <- grep("maxEAR",names(statsOfColumn))
  }
  

  freqCol <- grep("freq",names(statsOfColumn))
  n <- length(maxEARS)
  ignoreIndex <- which(names(statsOfColumn) %in% c("site","category"))
  
  if(n > 20 & n<30){
    colors <- c(brewer.pal(n = 12, name = "Set3"),
                  brewer.pal(n = 8, name = "Set2"),
                  brewer.pal(n = max(c(3,n-20)), name = "Set1"))
  } else if (n <= 20){
    colors <- c(brewer.pal(n = 12, name = "Set3"),
                  brewer.pal(n =  max(c(3,n-12)), name = "Set2"))     
  } else {
    colors <- colorRampPalette(brewer.pal(11,"Spectral"))(n)
  }
  
  tableSumm <- DT::datatable(statsOfColumn, extensions = 'Buttons',
                             rownames = FALSE,
                             options = list(#dom = 'ft',
                                            dom = 'Bfrtip',
                                            buttons =
                                              list('colvis'),
                                            scrollX = TRUE,
                                            # pageLength = nrow(statsOfColumn),
                                            order=list(list(colToSort,'desc'))))

  tableSumm <- formatRound(tableSumm, names(statsOfColumn)[-ignoreIndex], 2)

  for(i in 1:length(maxEARS)){
    tableSumm <- formatStyle(tableSumm,
                             names(statsOfColumn)[maxEARS[i]],
                             backgroundColor = colors[i])
    tableSumm <- formatStyle(tableSumm,
                             names(statsOfColumn)[freqCol[i]],
                             backgroundColor = colors[i])

    tableSumm <- formatStyle(tableSumm, names(statsOfColumn)[maxEARS[i]],
                             background = styleColorBar(range(statsOfColumn[,names(statsOfColumn)[maxEARS[i]]],na.rm = TRUE), 'goldenrod'),
                             backgroundSize = '100% 90%',
                             backgroundRepeat = 'no-repeat',
                             backgroundPosition = 'center' )
    tableSumm <- formatStyle(tableSumm, names(statsOfColumn)[freqCol[i]],
                             background = styleColorBar(range(statsOfColumn[,names(statsOfColumn)[freqCol[i]]],na.rm = TRUE), 'wheat'),
                             backgroundSize = '100% 90%',
                             backgroundRepeat = 'no-repeat',
                             backgroundPosition = 'center')

  }
  
  return(tableSumm)
}


#' @export
#' @rdname rank_sites_DT
rank_sites <- function(chemicalSummary, 
                           category, 
                           hit_threshold = 0.1, 
                           mean_logic = FALSE){
  
  sumEAR <- nHits <- n <- calc <- value <- choice_calc <- ".dplyr"
  chnm <- Class <- Bio_category <- site <- EAR <- maxEAR <- ".dplyr"
  mean_logic <- as.character(mean_logic)
  match.arg(mean_logic, c("mean","max","noSum","TRUE","FALSE"))
  
  if(mean_logic %in% c("TRUE","mean")){
    mean_logic <- TRUE
  }
  if(mean_logic %in% c("FALSE","max")){
    mean_logic <- FALSE
  }
  
  siteToFind <- unique(chemicalSummary$shortName)
  
  if(category == "Chemical"){
    chemicalSummary <- mutate(chemicalSummary, category = chnm)
  } else if (category == "Chemical Class"){
    chemicalSummary <- mutate(chemicalSummary, category = Class)
  } else {
    chemicalSummary <- mutate(chemicalSummary, category = Bio_category)
  }
  
  chemicalSummary <- select(chemicalSummary, -Class, -Bio_category, -chnm)
  
  if(length(siteToFind) == 1){
    chemicalSummary$site <- chemicalSummary$category
  } else {
    chemicalSummary$site <- chemicalSummary$shortName
  }
  
  if(mean_logic == "noSum"){
    statsOfColumn <- chemicalSummary %>%
      group_by(site, category) %>%
      summarise(maxEAR = max(EAR),
                freq = sum(maxEAR > 0)/n()) %>%
      data.frame()    
  } else {
    statsOfColumn <- chemicalSummary %>%
      group_by(site, date, category) %>%
      summarise(sumEAR = sum(EAR),
                nHits = sum(sumEAR > hit_threshold)) %>%
      group_by(site, category) %>%
      summarise(maxEAR = ifelse(mean_logic, mean(sumEAR), max(sumEAR)),
                freq = sum(nHits > 0)/n()) %>%
      data.frame()
  }

  if(!(length(siteToFind) == 1)){
    statsOfColumn <- statsOfColumn %>%
      gather(calc, value, -site, -category) %>%
      unite(choice_calc, category, calc, sep=" ") %>%
      spread(choice_calc, value)        
  }
  colToSort <- 2
  if("nSamples" %in% names(statsOfColumn)){
    colToSort <- 3
  }
  
  freqCol <- grep("freq",names(statsOfColumn))
  maxEARS <- grep("maxEAR",names(statsOfColumn))
  
  ignoreIndex <- which(names(statsOfColumn) %in% c("site","nSamples"))
  
  statsOfColumn <- statsOfColumn[,c(ignoreIndex,c(maxEARS,freqCol)[order(c(maxEARS,freqCol))])]
  
  maxEARS <- grep("maxEAR",names(statsOfColumn))
  
  MaxEARSordered <- order(apply(statsOfColumn[,maxEARS, drop = FALSE], 2, max),decreasing = TRUE)
  
  if(length(maxEARS) != 1){
    statsOfColumn <- statsOfColumn[,c(ignoreIndex,interl(maxEARS[MaxEARSordered],(maxEARS[MaxEARSordered]-1)))]
  }
  
  freqCol <- grep("freq",names(statsOfColumn))
  maxEARS <- grep("maxEAR",names(statsOfColumn))
  
  if(isTRUE(mean_logic)){
    names(statsOfColumn)[maxEARS] <- gsub("max","mean",names(statsOfColumn)[maxEARS])
  } 
  
  statsOfColumn <- statsOfColumn[order(statsOfColumn[[colToSort]], decreasing = TRUE),]
  
  if(length(siteToFind) == 1){
    names(statsOfColumn)[which(names(statsOfColumn) == "site")] <- "category"
  }
  return(statsOfColumn)
}


interl <- function (a,b) {
  n <- min(length(a),length(b))
  p1 <- as.vector(rbind(a[1:n],b[1:n]))
  p2 <- c(a[-(1:n)],b[-(1:n)])
  c(p1,p2)
}