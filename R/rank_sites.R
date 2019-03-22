#' Rank sites by EAR
#' 
#' The \code{rank_sites_DT} (DT option) and \code{rank_sites} (data frame option) functions 
#' create tables with one row per site. Columns represent the maximum or mean EAR 
#' (depending on the mean_logic argument) for each category ("Chemical Class", 
#' "Chemical", or "Biological") and the frequency of the maximum or mean EAR 
#' exceeding a user specified hit_threshold.
#' 
#' The tables show slightly different results for a single site. Rather than multiple 
#' columns for categories, there is now 1 row per category (since the site is known).
#' 
#' @param chemical_summary Data frame from \code{\link{get_chemical_summary}}.
#' @param mean_logic Logical.  \code{TRUE} displays the mean sample from each site,
#' \code{FALSE} displays the maximum sample from each site.
#' @param sum_logic Logical. \code{TRUE} sums the EARs in a specified grouping,
#' \code{FALSE} does not. \code{FALSE} may be better for traditional benchmarks as
#' opposed to ToxCast benchmarks.
#' @param category Character. Either "Biological", "Chemical Class", or "Chemical".
#' @param hit_threshold Numeric threshold defining a "hit".
#' @export
#' 
#' @return data frame with one row per site, and the mas or mean EAR and frequency of 
#' hits based on the category.
#' 
#' @rdname rank_sites_DT
#' @importFrom stats median
#' @examples
#' # This is the example workflow:
#' path_to_tox <-  system.file("extdata", package="toxEval")
#' file_name <- "OWC_data_fromSup.xlsx"
#'
#' full_path <- file.path(path_to_tox, file_name)
#' 
#' tox_list <- create_toxEval(full_path)
#' 
#' ACC <- get_ACC(tox_list$chem_info$CAS)
#' ACC <- remove_flags(ACC)
#' 
#' cleaned_ep <- clean_endPoint_info(end_point_info)
#' filtered_ep <- filter_groups(cleaned_ep)
#' chemical_summary <- get_chemical_summary(tox_list, ACC, filtered_ep)
#'
#' stats_df <- rank_sites(chemical_summary, "Biological")
#' 
#' rank_sites_DT(chemical_summary, category = "Biological")
#' rank_sites_DT(chemical_summary, category = "Chemical Class")
#' rank_sites_DT(chemical_summary, category = "Chemical")
#' 
rank_sites_DT <- function(chemical_summary, 
                          category = "Biological",
                          mean_logic = FALSE,
                          sum_logic = TRUE,
                          hit_threshold = 0.1){
  
  Bio_category <- Class <- EAR <- maxEAR <- sumEAR <- value <- calc <- chnm <- choice_calc <- n <- nHits <- site <- ".dplyr"
  
  match.arg(category, c("Biological","Chemical Class","Chemical"))

  statsOfColumn <- rank_sites(chemical_summary=chemical_summary,
                               category = category,
                               hit_threshold = hit_threshold,
                               mean_logic = mean_logic,
                               sum_logic = sum_logic)

  colToSort <- 1
  
  if(mean_logic){
    maxEARS <- grep("meanEAR",names(statsOfColumn))
  } else {
    maxEARS <- grep("maxEAR",names(statsOfColumn))
  }

  freqCol <- grep("freq",names(statsOfColumn))
  n <- length(maxEARS)
  ignoreIndex <- which(names(statsOfColumn) %in% c("site","category"))
  
  if(n > 20 & n<30){
    colors <- c(RColorBrewer::brewer.pal(n = 12, name = "Set3"),
                RColorBrewer::brewer.pal(n = 8, name = "Set2"),
                RColorBrewer::brewer.pal(n = max(c(3,n-20)), name = "Set1"))
  } else if (n <= 20){
    colors <- c(RColorBrewer::brewer.pal(n = 12, name = "Set3"),
                RColorBrewer::brewer.pal(n =  max(c(3,n-12)), name = "Set2"))     
  } else {
    colors <- colorRampPalette(RColorBrewer::brewer.pal(11,"Spectral"))(n)
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

  tableSumm <- DT::formatSignif(tableSumm, names(statsOfColumn)[-ignoreIndex], 2)

  for(i in 1:length(maxEARS)){
    tableSumm <- DT::formatStyle(tableSumm,
                             names(statsOfColumn)[maxEARS[i]],
                             backgroundColor = colors[i])
    tableSumm <- DT::formatStyle(tableSumm,
                             names(statsOfColumn)[freqCol[i]],
                             backgroundColor = colors[i])

    tableSumm <- DT::formatStyle(tableSumm, names(statsOfColumn)[maxEARS[i]],
                             background = DT::styleColorBar(range(statsOfColumn[,names(statsOfColumn)[maxEARS[i]]],na.rm = TRUE), 'goldenrod'),
                             backgroundSize = '100% 90%',
                             backgroundRepeat = 'no-repeat',
                             backgroundPosition = 'center' )
    tableSumm <- DT::formatStyle(tableSumm, names(statsOfColumn)[freqCol[i]],
                             background = DT::styleColorBar(range(statsOfColumn[,names(statsOfColumn)[freqCol[i]]],na.rm = TRUE), 'wheat'),
                             backgroundSize = '100% 90%',
                             backgroundRepeat = 'no-repeat',
                             backgroundPosition = 'center')

  }
  
  return(tableSumm)
}


#' @export
#' @rdname rank_sites_DT
rank_sites <- function(chemical_summary, 
                       category, 
                       hit_threshold = 0.1, 
                       mean_logic = FALSE,
                       sum_logic = TRUE){
  
  sumEAR <- nHits <- n <- calc <- value <- choice_calc <- ".dplyr"
  chnm <- Class <- Bio_category <- site <- EAR <- maxEAR <- ".dplyr"

  siteToFind <- unique(chemical_summary$shortName)
  
  if(category == "Chemical"){
    chemical_summary <- dplyr::mutate(chemical_summary, category = chnm)
  } else if (category == "Chemical Class"){
    chemical_summary <- dplyr::mutate(chemical_summary, category = Class)
  } else {
    chemical_summary <- dplyr::mutate(chemical_summary, category = Bio_category)
  }
  
  chemical_summary <- dplyr::select(chemical_summary, -Class, -Bio_category, -chnm)
  
  if(length(siteToFind) == 1){
    chemical_summary$site <- chemical_summary$category
  } else {
    chemical_summary$site <- chemical_summary$shortName
  }

  if(!sum_logic){
    statsOfColumn <- chemical_summary %>%
      dplyr::group_by(site, date, category) %>%
      dplyr::summarize(sumEAR = max(EAR),
                nHits = sum(sumEAR > hit_threshold)) %>%
      dplyr::group_by(site, category) %>%
      dplyr::summarise(maxEAR = ifelse(mean_logic, mean(sumEAR), max(sumEAR)),
                freq = sum(nHits > 0)/dplyr::n()) %>%
      data.frame()    
  } else {
    statsOfColumn <- chemical_summary %>%
      dplyr::group_by(site, date, category) %>%
      dplyr::summarise(sumEAR = sum(EAR),
                nHits = sum(sumEAR > hit_threshold)) %>%
      dplyr::group_by(site, category) %>%
      dplyr::summarise(maxEAR = ifelse(mean_logic, mean(sumEAR), max(sumEAR)),
                freq = sum(nHits > 0)/dplyr::n()) %>%
      data.frame()
  }

  if(!(length(siteToFind) == 1)){
    statsOfColumn <- statsOfColumn %>%
      tidyr::gather(calc, value, -site, -category) %>%
      tidyr::unite(choice_calc, category, calc, sep=" ") %>%
      tidyr::spread(choice_calc, value)        
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