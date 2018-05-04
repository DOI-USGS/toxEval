#' Biological hits per category
#' 
#' The \code{hits_by_groupings_DT} (DT option) and \code{hits_by_groupings} (data frame 
#' option) functions create tables with one row per category("Biological", 
#' "Chemical", or "Chemical Class"). The columns indicate the "Biological" 
#' groupings. The values in the table signify how many sites have exceeded 
#' the hit_threshold for that particular "Biological"/category combination. 
#' If the user chooses "Biological" as the category, it is a simple 2-column 
#' table of "Biological" groupings and number of sites (nSites).
#' 
#' The tables show slightly different results for a single site, showing the 
#' number of samples with hits (instead of number of sites).
#' 
#' @param chemicalSummary data frame from \code{get_chemical_summary}
#' @param mean_logic logical.  TRUE takes the mean sample of each site,
#' FALSE takes the maximum sample of each site.
#' @param sum_logic logical. TRUE sums the EARs in a specified grouping,
#' FALSE does not. FALSE may be better for traditional benchmarks as
#' opposed to ToxCast benchmarks.
#' @param category either "Biological", "Chemical Class", or "Chemical"
#' @param hit_threshold numeric threshold defining a "hit"
#' @export
#' @rdname hits_by_groupings_DT
#' @import DT
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
#' site_df <- hits_by_groupings(chemicalSummary, category = "Biological")
#' hits_by_groupings_DT(chemicalSummary, category = "Biological")
#' hits_by_groupings_DT(chemicalSummary, category = "Chemical Class")
#' hits_by_groupings_DT(chemicalSummary, category = "Chemical")
#' }
hits_by_groupings_DT <- function(chemicalSummary, 
                           category = "Biological",
                           mean_logic = FALSE,
                           sum_logic = TRUE,
                           hit_threshold = 0.1){
  
  match.arg(category, c("Biological","Chemical Class","Chemical"))

  tableData <- hits_by_groupings(chemicalSummary=chemicalSummary, 
                              category = category,
                              mean_logic = mean_logic,
                              sum_logic = sum_logic,
                              hit_threshold = hit_threshold)
  
  cuts <- seq(0,max(as.matrix(tableData),na.rm=TRUE),length.out = 8)
  colors <- brewer.pal(9,"Blues") #"RdYlBu"
  
  tableData1 <- DT::datatable(tableData, extensions = 'Buttons',
                              rownames = TRUE,
                              options = list(scrollX = TRUE,
                                             dom = 'Bfrtip',
                                             buttons = list('colvis'),
                               order=list(list(1,'desc'))))
  
  if(category != "Biological"){
    for(i in 1:ncol(tableData)){
      tableData1 <- formatStyle(tableData1, columns = names(tableData)[i],
                  backgroundColor = styleInterval(cuts = cuts,values = colors),
                  color = styleInterval(0.75*max(tableData,na.rm=TRUE),values = c("black","white")),
                  `font-size` = '17px')
    }
  }
  
  return(tableData1)
}

#' @export
#' @rdname hits_by_groupings_DT
hits_by_groupings <- function(chemicalSummary, 
                              category, 
                              mean_logic=FALSE, 
                              sum_logic=TRUE,
                              hit_threshold = 0.1){
  
  Bio_category <- Class <- EAR <- sumEAR <- value <- calc <- chnm <- choice_calc <- n <- nHits <- site <- ".dplyr"
  meanEAR <- nSites <- ".dplyr"
  match.arg(category, c("Biological","Chemical Class","Chemical"))

  if(category == "Biological"){
    chemicalSummary$category <- chemicalSummary$Bio_category
  } else if(category == "Chemical Class") {
    chemicalSummary$category <- chemicalSummary$Class
  } else {
    chemicalSummary$category <- chemicalSummary$chnm
  }
  
  if(length(unique(chemicalSummary$site)) > 1){
    
    if(!sum_logic){
      tableData <- chemicalSummary %>%
        group_by(site, Bio_category, category) %>%
        summarize(meanEAR = ifelse(mean_logic, mean(EAR),max(EAR))) %>%
        group_by(Bio_category, category) %>%
        summarize(nSites = sum(meanEAR >  hit_threshold)) %>%
        data.frame()
    } else {
      tableData <- chemicalSummary %>%
        group_by(site, Bio_category, category, date) %>%
        summarize(sumEAR = sum(EAR)) %>%
        group_by(site, Bio_category, category) %>%
        summarize(meanEAR = ifelse(mean_logic, mean(sumEAR),max(sumEAR))) %>%
        group_by(Bio_category, category) %>%
        summarize(nSites = sum(meanEAR >  hit_threshold)) %>%
        data.frame()      
    }
  } else {
    
    if(!sum_logic){
      tableData <- chemicalSummary %>%
        group_by(Bio_category, category) %>%
        summarise(nSites = sum(EAR > hit_threshold))%>%
        data.frame()      
    } else {
      tableData <- chemicalSummary %>%
        group_by(Bio_category, category, date)%>%
        summarise(sumEAR=sum(EAR)) %>%
        data.frame() %>%
        group_by(Bio_category, category) %>%
        summarise(nSites = sum(sumEAR > hit_threshold))%>%
        data.frame()      
    }
  }
  
  if(category != "Biological"){
    
    tableData <- tableData %>%
      spread(Bio_category, nSites)
    
    sumOfColumns <- colSums(tableData[-1],na.rm = TRUE)
    
    if(!all(sumOfColumns == 0)){
      orderData <- order(sumOfColumns,decreasing = TRUE)
      orderData <- orderData[sumOfColumns[orderData] != 0] + 1
      
      tableData <- tableData[,c(1,orderData)]      
    }

    groups <- tableData$category
    
    tableData <- tableData[!is.na(groups),-1,drop=FALSE]
    rownames(tableData) <- groups[!is.na(groups)]
    
  } else {
    tableData <- select(tableData, Bio_category, nSites)
    rownames(tableData) <- tableData$Bio_category
    tableData <- tableData[,-1,drop=FALSE]
  }
  
  if(length(unique(chemicalSummary$site)) == 1){
    names(tableData)[names(tableData) == "nSites"] <- "nSamples" 
  }
  
  return(tableData)
}
