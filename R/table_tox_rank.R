#' table_tox_rank
#' 
#' Table of ranks
#' @param chemicalSummary data frame from \code{get_chemical_summary}
#' @param mean_logic logical \code{TRUE} is mean, \code{FALSE} is maximum
#' @param category either "Biological", "Chemical Class", or "Chemical"
#' @param hit_threshold numeric threshold defining a "hit"
#' @export
#' @import DT
#' @importFrom stats median
#' @importFrom tidyr spread unite
#' @importFrom dplyr full_join filter mutate select left_join right_join
#' @examples
#' library(readxl)
#' path_to_tox <-  system.file("extdata", package="toxEval")
#' file_name <- "OWC_data_fromSup.xlsx"
#' full_path <- file.path(path_to_tox, file_name)
#' 
#' chem_data <- read_excel(full_path, sheet = "Data")
#' chem_info <- read_excel(full_path, sheet = "Chemicals") 
#' chem_site <- read_excel(full_path, sheet = "Sites")
#' ACClong <- get_ACC(chem_info$CAS)
#' ACClong <- remove_flags(ACClong)
#' 
#' cleaned_ep <- clean_endPoint_info(endPointInfo)
#' filtered_ep <- filter_groups(cleaned_ep)
#' 
#' chemicalSummary <- get_chemical_summary(ACClong,
#'                                         filtered_ep,
#'                                        chem_data, 
#'                                         chem_site, 
#'                                         chem_info)
#' table_tox_rank(chemicalSummary, category = "Biological")
#' table_tox_rank(chemicalSummary, category = "Chemical Class")
#' table_tox_rank(chemicalSummary, category = "Chemical")
table_tox_rank <- function(chemicalSummary, 
                          category = "Biological",
                          mean_logic = FALSE,
                          hit_threshold = 0.1){
  
  Bio_category <- Class <- EAR <- sumEAR <- value <- calc <- chnm <- choice_calc <- n <- nHits <- site <- ".dplyr"
  
  match.arg(category, c("Biological","Chemical Class","Chemical"))
  
  siteToFind <- unique(chemicalSummary$shortName)

  statsOfColumn <- statsOfColumns(chemicalSummary=chemicalSummary,
                                   category = category,
                                   hit_threshold = hit_threshold,
                                   mean_logic = mean_logic)

  colToSort <- 1
  if("nSamples" %in% names(statsOfColumn)){
    colToSort <- 2
  }
  freqCol <- grep("freq",names(statsOfColumn))
  maxEARS <- grep("maxEAR",names(statsOfColumn))

  ignoreIndex <- which(names(statsOfColumn) %in% c("site","nSamples"))

  statsOfColumn <- statsOfColumn[,c(ignoreIndex,c(maxEARS,freqCol)[order(c(maxEARS,freqCol))])]

  maxEARS <- grep("maxEAR",names(statsOfColumn))

  MaxEARSordered <- order(apply(statsOfColumn[,maxEARS, drop = FALSE], 2, max),decreasing = TRUE)

  if(length(maxEARS) > 9){
    statsOfColumn <- statsOfColumn[,c(ignoreIndex,interl(maxEARS[MaxEARSordered[1:9]],(maxEARS[MaxEARSordered[1:9]]-1)))]
    maxEARS <- maxEARS[1:9]
  } else {
    if(length(maxEARS) != 1){
      statsOfColumn <- statsOfColumn[,c(ignoreIndex,interl(maxEARS[MaxEARSordered],(maxEARS[MaxEARSordered]-1)))]
    }
  }

  freqCol <- grep("freq",names(statsOfColumn))
  maxEARS <- grep("maxEAR",names(statsOfColumn))

  if(mean_logic){
    names(statsOfColumn)[maxEARS] <- gsub("max","mean",names(statsOfColumn)[maxEARS])
  }

  colors <- brewer.pal(length(maxEARS),"Blues") #"RdYlBu"
  
  tableSumm <- DT::datatable(statsOfColumn, extensions = 'Buttons',
                             rownames = FALSE,
                             options = list(#dom = 'ft',
                                            dom = 'Bfrtip',
                                            buttons =
                                              list('colvis', list(
                                                extend = 'collection',
                                                buttons = list(list(extend='csv',
                                                                    filename = 'hitStats'),
                                                               list(extend='excel',
                                                                    filename = 'hitStats'),
                                                               list(extend='pdf',
                                                                    filename= 'hitStats')),
                                                text = 'Download'
                                                )
                                              ),
                                            scrollX = TRUE,
                                            pageLength = nrow(statsOfColumn),
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

#' statsOfColumns
#' 
#' Summarize data for most graphs/tables
#' @param chemicalSummary data frame
#' @param category character
#' @param hit_threshold numeric
#' @param mean_logic logical
#' @export
#' @keywords internal
statsOfColumns <- function(chemicalSummary, 
                           category, 
                           hit_threshold, 
                           mean_logic){
  
  sumEAR <- nHits <- n <- calc <- value <- choice_calc <- ".dplyr"
  
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
  
  statsOfColumn <- chemicalSummary %>%
    group_by(site, date, category) %>%
    summarise(sumEAR = sum(EAR),
              nHits = sum(EAR > hit_threshold)) %>%
    group_by(site, category) %>%
    summarise(maxEAR = ifelse(mean_logic, mean(sumEAR), max(sumEAR)),
              freq = sum(nHits > 0)/n()) %>%
    data.frame()
  
  if(!(length(siteToFind) == 1)){
    statsOfColumn <- statsOfColumn %>%
      gather(calc, value, -site, -category) %>%
      unite(choice_calc, category, calc, sep=" ") %>%
      spread(choice_calc, value)        
  }
  
  return(statsOfColumn)
}

#' interl
#' 
#' Interleaves vector
#' @export
#' @keywords internal
#' @param a vector
#' @param b vector
interl <- function (a,b) {
  n <- min(length(a),length(b))
  p1 <- as.vector(rbind(a[1:n],b[1:n]))
  p2 <- c(a[-(1:n)],b[-(1:n)])
  c(p1,p2)
}