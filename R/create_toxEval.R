#' create_toxEval
#' 
#' 
#' @param excel_file_path Path to Excel file that contains at least 3 tabs: Data, Chemicals, and Sites, 
#' and could optionally contain Exclude and Benchmark
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
    stop("Excel file needs at least 3 tabs labeled 'Data', 'Chemicals', and 'Sites'")
  }
  
  req_cols <- c("CAS", "SiteID", "Value", "Sample Date")
  chem_data <- get_and_check("Data",req_cols, excel_file_path)

  req_cols <- c("CAS", "Class")
  chem_info <- get_and_check("Chemicals",req_cols, excel_file_path)

  req_cols <- c("SiteID", "Short Name")
  chem_site <- get_and_check("Sites",req_cols, excel_file_path)
  
  exclusions <- NULL
  benchmarks <- NULL
  
  if("Exclude" %in% readxl::excel_sheets(excel_file_path)){
    req_cols <- c("CAS", "endPoint")
    exclusions <- get_and_check("Exclude",req_cols, excel_file_path)
  }
  
  if("Benchmarks" %in% readxl::excel_sheets(excel_file_path)){
    req_cols <- c("CAS", "endPoint")
    benchmarks <- get_and_check("Benchmarks",req_cols, excel_file_path)
  }
  
  #Check that all CAS in Data in Chemical
  data_cas <- unique(chem_data$CAS)
  chem_cas <- unique(chem_info$CAS)
  if(!all(data_cas %in% chem_cas)){
    warning("Not all CAS in the Data tab are listed in the Chemical tab")
  }
  
  #Check that all SiteID in Data in Sites
  data_siteID <- unique(chem_data$SiteID)
  site_siteID <- unique(chem_site$SiteID)
  if(!all(data_siteID %in% site_siteID)){
    warning("Not all SiteID in the Data tab are listed in the Sites tab")
  }
  rawData <- list(chem_data=chem_data,
                  chem_info=chem_info,
                  chem_site=chem_site,
                  exclusions=exclusions,
                  benchmarks=benchmarks)
  
  return(rawData)
  
}

get_and_check <- function(tab_name, req_cols, excel_file_path){

  tab <- readxl::read_excel(excel_file_path, sheet = tab_name)
  if(!all(req_cols %in% names(tab))){
    miss_col <- req_cols[!req_cols %in% names(tab)]
    stop(tab," tab missing ",miss_col," column")
  }
  
  return(tab)
}