#' Load and check toxEval data
#' 
#' This function is used to load a data file for analysis in the form of a 
#' single Excel file. The Excel file should include 3 mandatory sheets named 
#' "Data", "Chemicals", and "Sites". Additionally there are 2 optional sheets: 
#' "Exclude" and "Benchmarks". This function creates a data frame for each 
#' sheet,  perform basic checks on the data to assure that required columns are
#' included for each sheet 
#' 
#' Required columns in the Data sheet include "CAS", "SiteID", "Value", and 
#' "Sample Date". The "Value" column includes concentration measurements in 
#' units of \eqn{\mu}g/L. "Sample Date" can be either a date or date/time or an integer.
#' Additional columns may be included for user purposes, but will not be used in
#' toxEval. 
#' 
#' Required columns in the Chemical sheet include "CAS", "Class". "CAS" values in 
#' this sheet must exactly match corresponding "CAS" values in the Data sheet. The 
#' "Class" designation allows  data to be grouped in a user-specified way. For 
#' example, in a data set of multiple pesticides, it may be valuable to explore 
#' differences and similarities to of insecticides, herbicides and fungicides. 
#' Additional columns may be included for user purposes, but will not be used in
#' toxEval. 
#' 
#' Required columns in the Sites sheet includes "SiteID", "Short Name", and for 
#' the Shiny application "dec_lat","dec_lon". Values in the "SiteID" column in 
#' this sheet exactly match corresponding values in the "SiteID" column in the 
#' Data sheet. Additional columns may be included for user purposes, but will not 
#' be used in toxEval. 
#' 
#' When using the optional sheet Exclude, columns required include "CAS" and 
#' "endPoint". These are used to exclude particular chemicals (via CAS), 
#' ToxCast endpoints (via endPoint), or a unique chemical/endpoint combination. 
#' Additional columns may be included for user purposes, but will not be used 
#' in toxEval. 
#' 
#' When using the optional sheet Benchmarks, columns required include "CAS", 
#' "endPoint","ACC_value" and "chnm". This sheet is used to over-ride the 
#' functions using endpoints from the ToxCast database, allowing the user
#' to import endpoint information from other sources. It could also be useful 
#' for reproducing results in the future (for example, if after ToxCast updates, 
#' analysis with an older version of ToxCast may be reproduced by including the
#' selected ToxCast endpoint database in this sheet. Additional columns may be included for user 
#' purposes, but will not be used in toxEval. 
#' 
#' 
#' For more information, see the "Prepare Data" vignette: \url{../doc/PrepareData.html}.
#' 
#' All remaining toxEval functions use data from via the list that is returned 
#' from this function.
#'  
#' @return The object returned from this function contains a list of between 
#' three and five data frames. The minimum data frames returned are chem_data 
#' (containing at least the columns: "CAS", "SiteID", "Value", "Sample Date"),
#' chem_info (containing at least the columns: "CAS", "Class"), and
#' chem_site (containing at least the columns: "SiteID", "Short Name". For the
#' Shiny app, "dec_lat" and "dec_lon" are also required). The optional data 
#' frames are exclusions (containing at least the columns: "CAS" and "endPoint"),
#' and benchmarks (containing at least the columns: "CAS", "endPoint",
#' "ACC_value" and "chnm")
#' 
#' @param excel_file_path Path to Excel file that contains at least 3 sheets: Data, Chemicals, and Sites, 
#' and could optionally contain Exclude and Benchmarks.
#'
#' @export
#' @examples 
#' path_to_tox <-  system.file("extdata", package="toxEval")
#' file_name <- "OWC_data_fromSup.xlsx"
#' excel_file_path <- file.path(path_to_tox, file_name)
#' tox_list <- create_toxEval(excel_file_path)
create_toxEval <- function(excel_file_path){
  
  if(!file.exists(excel_file_path)){
    stop("File does not exist, check path and spelling")
  } 
  
  if(!all(c("Data","Chemicals","Sites") %in% readxl::excel_sheets(excel_file_path))){
    stop("Excel file needs at least 3 sheets labeled 'Data', 'Chemicals', and 'Sites'")
  }
  
  req_cols <- c("CAS", "SiteID", "Value", "Sample Date")
  chem_data <- readxl::read_excel(excel_file_path, sheet = "Data")
  # Correct some common possibilities:
  names(chem_data)[names(chem_data) %in% c("site","Site","Station","station")] <- "SiteID"
  names(chem_data)[names(chem_data) %in% c("Date","date","sample","Sample")] <- "Sample Date"
  names(chem_data)[names(chem_data) %in% c("casrn", "casn","CASRN","CASN")] <- "CAS"
  chem_data <- check_cols(req_cols, chem_data)
  
  if(!is.numeric(chem_data$Value)){
    stop("The 'Value' column in the Data sheet is not numeric")
  }
  
  req_cols <- c("CAS", "Class")
  chem_info <- readxl::read_excel(excel_file_path, sheet = "Chemicals")
  names(chem_info)[names(chem_info) %in% c("casrn", "casn","CASRN","CASN")] <- "CAS"
  chem_info <- check_cols(req_cols, chem_info)

  req_cols <- c("SiteID", "Short Name")
  chem_site <- readxl::read_excel(excel_file_path, sheet = "Sites")
  names(chem_site)[names(chem_site) %in% c("site","Site","Station","station")] <- "SiteID"
  names(chem_site)[names(chem_site) %in% c("shortName", "map_name")] <- "Short Name" #others?
  names(chem_site)[names(chem_site) %in% c("dec_long", "longitude", "long")] <- "dec_lon"
  names(chem_site)[names(chem_site) %in% c("lat", "latitude")] <- "dec_lat"
  chem_site <- check_cols(req_cols, chem_site)
  
  #Check columns needed for shiny app:
  app_cols <- c("dec_lat","dec_lon")
  chem_site <- check_cols(app_cols, chem_site, mandatory = FALSE)
  
  if(!("site_grouping" %in% names(chem_site))){
    chem_site$site_grouping <- ""
  }
  
  exclusions <- NULL
  benchmarks <- NULL
  
  if("Exclude" %in% readxl::excel_sheets(excel_file_path)){
    req_cols <- c("CAS", "endPoint")
    exclusions <- readxl::read_excel(excel_file_path, sheet = "Exclude")
    names(exclusions)[names(exclusions) %in% c("casrn", "casn","CASRN","CASN")] <- "CAS"
    exclusions <- check_cols(req_cols, exclusions)
  }
  
  if("Benchmarks" %in% readxl::excel_sheets(excel_file_path)){
    req_cols <- c("CAS", "endPoint","ACC_value","chnm")
    benchmarks <- readxl::read_excel(excel_file_path, sheet = "Benchmarks")
    
    names(benchmarks)[names(benchmarks) %in% c("casrn", "casn","CASRN","CASN")] <- "CAS"
    names(benchmarks)[names(benchmarks) %in% c("Value", "value","ACC","ACC_value")] <- "ACC_value"
    names(benchmarks)[names(benchmarks) %in% c("chemical", "chnm","Chemical","Compound")] <- "chnm"
    
    benchmarks <- check_cols(req_cols, benchmarks)
    if(!("groupCol" %in% names(benchmarks))){
      benchmarks$groupCol <- "Benchmarks"
    }
  }
  
  #Check that all CAS in Data in Chemical
  data_cas <- unique(chem_data$CAS)
  chem_cas <- unique(chem_info$CAS)
  if(!all(data_cas %in% chem_cas)){
    warning("Not all CAS in the Data sheet are listed in the Chemical sheet")
  }
  
  #Check that all SiteID in Data in Sites
  data_siteID <- unique(chem_data$SiteID)
  site_siteID <- unique(chem_site$SiteID)
  if(!all(data_siteID %in% site_siteID)){
    warning("Not all SiteID in the Data sheet are listed in the Sites sheet")
  }
  rawData <- list(chem_data=chem_data,
                  chem_info=chem_info,
                  chem_site=chem_site,
                  exclusions=exclusions,
                  benchmarks=benchmarks)
  
  class(rawData) <- "toxEval"
  
  return(rawData)
  
}

check_cols <- function(req_cols, tab, mandatory=TRUE){
  
  #allow column names to be case and whitespace insensitive:
  req_cols_low <- tolower(req_cols)
  req_cols_low <- gsub(" ", "", req_cols_low)
  req_cols_low <- gsub("_", "", req_cols_low)
  req_cols_low <- gsub("\\.", "", req_cols_low)
  
  orig_names <- names(tab)
  names(tab) <- tolower(names(tab))
  names(tab) <- gsub(" ", "", names(tab))
  names(tab) <- gsub("_", "", names(tab))
  names(tab) <- gsub("\\.", "", names(tab))
  
  if(!all(req_cols_low %in% names(tab))){
    miss_col <- req_cols_low[!req_cols_low %in% names(tab)]
    
    if(mandatory){
      stop(tab," tab missing ",miss_col," column")
    } else {
      warning(tab," tab missing optional column: ",miss_col,".\n
              This could cause Shiny app to not work properly")
    }
  }
  # Transform the column names the "required" names:
  names(tab)[!(names(tab) %in% req_cols_low)] <- orig_names[-which(names(tab) %in% req_cols_low)]
  names(tab)[(names(tab) %in% req_cols_low)] <- unlist(sapply(names(tab), function(x) req_cols[which(req_cols_low == x)]))
  
  tab <- rm_em_dash(tab)

  return(tab)
}


rm_em_dash <- function(df){
  for(i in as.integer(which(sapply(df, class) == "character"))){
    df[[i]] <- gsub("\u002D|\u2013|\u2014|\u2212","-",df[[i]])
  }
  return(df)
}


#' summary of tox_list
#' 
#' 
#' @param object toxEval object with "chem_info" data frame included. 
#' @param \dots additional parameters
#'
#' @export
#' @examples 
#' path_to_tox <-  system.file("extdata", package="toxEval")
#' file_name <- "OWC_data_fromSup.xlsx"
#' excel_file_path <- file.path(path_to_tox, file_name)
#' tox_list <- create_toxEval(excel_file_path)
#' summary(tox_list) 
summary.toxEval <- function(object, ...){
  
  casn <- endPoint <- chnm <- flags <- ".dplyr"
  
  if(is.null(object[["benchmarks"]])){
    ACClong <- ACC %>%
      dplyr::filter(casn %in% unique(object$chem_info$CAS)) %>%
      tidyr::gather(endPoint, ACC, -casn, -chnm, -flags) 
    bench_word <- "ToxCast"
  } else {
    ACClong <- object[["benchmarks"]] 
    bench_word <- "benchmark"
  }
 
  CAS_tots <- dplyr::select(ACClong, casn) %>% dplyr::distinct() %>% dplyr::pull(casn)
  
  CAS_w_data <- ACClong %>% dplyr::select(ACC, casn) %>%
    dplyr::filter(!is.na(ACC)) %>%
    dplyr::select(casn) %>%
    dplyr::distinct() %>% dplyr::pull(casn)
  
  message(length(CAS_tots)," chemicals have ", bench_word, " information")
  message("Chemicals returned from this function do NOT have ", bench_word, " information:")
  return(unique(object$chem_info$CAS)[!(unique(object$chem_info$CAS) %in% CAS_tots)])
}