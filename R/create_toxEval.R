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
#' For more information, see the "Prepare Data" vignette: \href{../doc/PrepareData.html}{\code{vignette("PrepareData", package = "toxEval")}}.
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
#' @param \dots data frames to override data within the original x list.
#'
#' @export
#' @examples
#'
#' path_to_tox <- system.file("extdata", package = "toxEval")
#' file_name <- "OWC_data_fromSup.xlsx"
#' excel_file_path <- file.path(path_to_tox, file_name)
#' tox_list <- create_toxEval(excel_file_path)
create_toxEval <- function(excel_file_path, ...) {
  if (!file.exists(excel_file_path)) {
    stop("File does not exist, check path and spelling")
  }

  if ("Chemicals" %in% readxl::excel_sheets(excel_file_path)) {
    chem_info <- readxl::read_excel(excel_file_path,
      sheet = "Chemicals",
      guess_max = 10000
    )

    names(chem_info) <- names(chem_info) %>%
      allowed_names(c("casrn", "casn", "CASRN", "CASN"), "CAS")

    if (all(c("CAS", "Chemical") %in% names(chem_info))){
      chem_info$CAS <- as.character(chem_info$CAS)
    } else {
      message("Chemical tab missing CAS or Chemical column")
    }
  }

  if ("Data" %in% readxl::excel_sheets(excel_file_path)) {
    chem_data <- readxl::read_excel(excel_file_path,
      sheet = "Data",
      guess_max = 100000
    )

    names(chem_data) <- names(chem_data) %>%
      allowed_names(c("site", "Site", "Station", "station"), "SiteID") %>%
      allowed_names(c("Date", "date", "sample", "Sample"), "Sample Date") %>%
      allowed_names(c("casrn", "casn", "CASRN", "CASN"), "CAS")

    if ("CAS" %in% names(chem_data)) {
      chem_data$CAS <- as.character(chem_data$CAS)
    } else {
      message("Data tab missing CAS column")
    }

    if ("SiteID" %in% names(chem_data)) {
      chem_data$SiteID <- as.character(chem_data$SiteID)
    } else {
      message("Data tab missing SiteID column")
    }
  }

  if ("Sites" %in% readxl::excel_sheets(excel_file_path)) {
    chem_site <- readxl::read_excel(excel_file_path,
      sheet = "Sites",
      guess_max = 5000
    )

    names(chem_site) <- names(chem_site) %>%
      allowed_names(c("site", "Site", "Station", "station"), "SiteID") %>%
      allowed_names(c("shortName", "map_name"), "Short Name") %>%
      allowed_names(c("dec_long", "longitude", "long"), "dec_lon") %>%
      allowed_names(c("lat", "latitude"), "dec_lat")

    if ("SiteID" %in% names(chem_site)) {
      chem_site$SiteID <- as.character(chem_site$SiteID)
    } else {
      message("Sites tab missing SiteID column")
    }
  }

  if ("Exclude" %in% readxl::excel_sheets(excel_file_path)) {
    exclusions <- readxl::read_excel(excel_file_path,
      sheet = "Exclude",
      guess_max = 1000
    )
  } else {
    exclusions <- NULL
  }

  if ("Benchmarks" %in% readxl::excel_sheets(excel_file_path)) {
    benchmarks <- readxl::read_excel(excel_file_path,
      sheet = "Benchmarks",
      guess_max = 10000
    )
    if (nrow(benchmarks) == 0) {
      benchmarks <- NULL
    }
  } else {
    benchmarks <- NULL
  }

  rawTox <- list(
    chem_data = chem_data,
    chem_info = chem_info,
    chem_site = chem_site,
    exclusions = exclusions,
    benchmarks = benchmarks
  )

  rawTox <- as.toxEval(rawTox, ...)

  return(rawTox)
}

check_cols <- function(req_cols, tab, tab_name, mandatory = TRUE) {

  # allow column names to be case and whitespace insensitive:
  req_cols_low <- tolower(req_cols)
  req_cols_low <- gsub(" ", "", req_cols_low)
  req_cols_low <- gsub("_", "", req_cols_low)
  req_cols_low <- gsub("\\.", "", req_cols_low)

  orig_names <- names(tab)
  names(tab) <- tolower(names(tab))
  names(tab) <- gsub(" ", "", names(tab))
  names(tab) <- gsub("_", "", names(tab))
  names(tab) <- gsub("\\.", "", names(tab))

  if (!all(req_cols_low %in% names(tab))) {
    miss_col <- req_cols_low[!req_cols_low %in% names(tab)]

    if (mandatory) {
      stop("Missing ", miss_col, " column in ", tab_name)
    } else {
      warning(
        "Missing optional column: ", miss_col, " in", tab_name,
        ".\nThis could cause Shiny app to not work properly"
      )
    }
  }
  # Transform the column names the "required" names:
  names(tab)[!(names(tab) %in% req_cols_low)] <- orig_names[-which(names(tab) %in% req_cols_low)]
  names(tab)[(names(tab) %in% req_cols_low)] <- unlist(sapply(names(tab), function(x) req_cols[which(req_cols_low == x)]))

  tab <- rm_em_dash(tab)

  return(tab)
}


rm_em_dash <- function(df) {
  for (i in as.integer(which(sapply(df, class) == "character"))) {
    df[[i]] <- gsub("\u002D|\u2013|\u2014|\u2212", "-", df[[i]])
  }
  return(df)
}


#' Summary of tox_list
#'
#' A "tox_list" object is created from \code{create_toxEval}. It
#' is a list of  5 data frames: chem_data, chem_info,
#' chem_site, exclusions, and benchmarks. This function returns
#' a message with how many chemicals have ToxCast information,
#' and returns a vector of which chemicals do not have ToxCast information.
#'
#' @param object toxEval object with "chem_info" data frame included.
#' @param \dots additional parameters
#'
#' @export
#' @importFrom magrittr "%>%"
#' @examples
#' path_to_tox <- system.file("extdata", package = "toxEval")
#' file_name <- "OWC_data_fromSup.xlsx"
#' excel_file_path <- file.path(path_to_tox, file_name)
#' tox_list <- create_toxEval(excel_file_path)
#' summary(tox_list)
summary.toxEval <- function(object, ...) {
  CAS <- endPoint <- chnm <- flags <- ".dplyr"

  if (is.null(object[["benchmarks"]])) {
    ACC <- ToxCast_ACC %>%
      filter(CAS %in% unique(object$chem_info$CAS))
    bench_word <- "ToxCast"
  } else {
    ACC <- object[["benchmarks"]]
    bench_word <- "benchmark"
  }

  CAS_w_data <- ACC %>%
    select(CAS) %>%
    distinct() %>%
    pull(CAS)

  message(length(CAS_w_data), " chemicals have ", bench_word, " information")
  message("Chemicals returned from this function do NOT have ", bench_word, " information:")
  return(unique(object$chem_info$CAS)[!(unique(object$chem_info$CAS) %in% CAS_w_data)])
}


allowed_names <- function(col_names, allowable, better) {
  if (any(allowable %in% col_names)) {
    col_names[col_names %in% allowable] <- better
  }
  return(col_names)
}


#' toxEval helper functions
#'
#' A small collection of helper functions for toxEval
#'
#' @export
#' @param x list or toxEval object
#' @param \dots data frames to override data within the original x list.
#' @rdname helperToxEval
#'
#' @examples
#' path_to_tox <- system.file("extdata", package = "toxEval")
#' file_name <- "OWC_data_fromSup.xlsx"
#' full_path <- file.path(path_to_tox, file_name)
#' tox_list <- create_toxEval(full_path)
#'
#' # To over-ride one of the data frames:
#' chem_data <- data.frame(
#'   CAS = "21145-77-7",
#'   Value = 1,
#'   "Sample Date" = as.Date("2012-01-01"),
#'   SiteID = "USGS-04249000"
#' )
#' tox_list_new <- as.toxEval(tox_list, chem_data)
#'
as.toxEval <- function(x, ...) {
  matchReturn <- c(
    do.call("c", list(...)[sapply(list(...), class) == "list"]), # get the list parts
    list(...)[sapply(list(...), class) != "list"]
  )

  if (is.null(names(matchReturn))) {
    names(matchReturn) <- match.call(expand.dots = FALSE)$...
  }

  if ("chem_data" %in% names(matchReturn)) {
    message("Replacing chem_data from list with supplied data frame.")
    x[["chem_data"]] <- matchReturn[["chem_data"]]
  }

  if ("chem_info" %in% names(matchReturn)) {
    message("Replacing chem_info from list with supplied data frame.")
    x[["chem_info"]] <- matchReturn[["chem_info"]]
  }

  if ("chem_site" %in% names(matchReturn)) {
    message("Replacing chem_site from list with supplied data frame.")
    x[["chem_site"]] <- matchReturn[["chem_site"]]
  }

  if ("exclusions" %in% names(matchReturn)) {
    message("Replacing exclusions from list with supplied data frame.")
    x[["exclusions"]] <- matchReturn[["exclusions"]]
  }

  if ("benchmarks" %in% names(matchReturn)) {
    message("Replacing benchmarks from list with supplied data frame.")
    x[["benchmarks"]] <- matchReturn[["benchmarks"]]
  }

  names(x) <- names(x) %>%
    tolower() %>%
    allowed_names("data", "chem_data") %>%
    allowed_names("chemicals", "chem_info") %>%
    allowed_names("sites", "chem_site") %>%
    allowed_names("exclude", "exclusions")

  required_data <- c("chem_data", "chem_info", "chem_site")

  if (!all(required_data %in% names(x))) {
    stop("Missing", required_data[!required_data %in% names(x)])
  }

  chem_data <- x[["chem_data"]]
  chem_info <- x[["chem_info"]]
  chem_site <- x[["chem_site"]]

  if ("benchmarks" %in% names(x)) {
    benchmarks <- x[["benchmarks"]]
  } else {
    benchmarks <- NULL
  }

  if ("exclusions" %in% names(x)) {
    exclusions <- x[["exclusions"]]
  } else {
    exclusions <- NULL
  }

  # chem_data:
  names(chem_data) <- names(chem_data) %>%
    allowed_names(c("site", "Site", "Station", "station"), "SiteID") %>%
    allowed_names(c("Date", "date", "sample", "Sample"), "Sample Date") %>%
    allowed_names(c("casrn", "casn", "CASRN", "CASN"), "CAS")

  chem_data <- check_cols(c("CAS", "SiteID", "Value", "Sample Date"), chem_data, tab_name = "chem_data")

  if (!is.numeric(chem_data$Value)) {
    stop("The 'Value' column in the Data sheet is not numeric")
  }

  # chem_info:
  names(chem_info) <- names(chem_info) %>%
    allowed_names(c("casrn", "casn", "CASRN", "CASN"), "CAS")

  chem_info <- check_cols(c("CAS", "Class"), chem_info, tab_name = "chem_info")

  # chem_site:
  names(chem_site) <- names(chem_site) %>%
    allowed_names(c("site", "Site", "Station", "station"), "SiteID") %>%
    allowed_names(c("shortName", "map_name"), "Short Name") %>%
    allowed_names(c("dec_long", "longitude", "long"), "dec_lon") %>%
    allowed_names(c("lat", "latitude"), "dec_lat")

  chem_site <- check_cols(c("SiteID", "Short Name"), chem_site, tab_name = "chem_site")

  # Check columns needed for shiny app:
  chem_site <- check_cols(c("dec_lat", "dec_lon"),
    chem_site,
    mandatory = FALSE, tab_name = "chem_site"
  )

  if (!("site_grouping" %in% names(chem_site))) {
    chem_site$site_grouping <- ""
  }

  if (!is.null(exclusions)) {
    names(exclusions) <- names(exclusions) %>%
      allowed_names(c("casrn", "casn", "CASRN", "CASN"), "CAS")

    exclusions <- check_cols(c("CAS", "endPoint"), exclusions, tab_name = "exclusions")
  }

  if (!is.null(benchmarks)) {
    names(benchmarks) <- names(benchmarks) %>%
      allowed_names(c("casrn", "casn", "CASRN", "CASN"), "CAS") %>%
      allowed_names(c("Value", "value", "ACC", "ACC_value"), "ACC_value") %>%
      allowed_names(c("chemical", "chnm", "Chemical", "Compound"), "chnm")

    benchmarks <- check_cols(c("CAS", "endPoint", "ACC_value", "chnm"),
      benchmarks,
      tab_name = "benchmarks"
    )
    if (!("groupCol" %in% names(benchmarks))) {
      benchmarks$groupCol <- "Benchmarks"
    }
  }

  # Check that all CAS in Data in Chemical
  data_cas <- unique(chem_data$CAS)
  chem_cas <- unique(chem_info$CAS)
  if (!all(data_cas %in% chem_cas)) {
    warning("Not all CAS in the Data sheet are listed in the Chemical sheet")
  }

  # Check that all SiteID in Data in Sites
  data_siteID <- unique(chem_data$SiteID)
  site_siteID <- unique(chem_site$SiteID)
  if (!all(data_siteID %in% site_siteID)) {
    warning("Not all SiteID in the Data sheet are listed in the Sites sheet")
  }
  rawData <- list(
    chem_data = chem_data,
    chem_info = chem_info,
    chem_site = chem_site,
    exclusions = exclusions,
    benchmarks = benchmarks
  )

  class(rawData) <- "toxEval"

  return(rawData)
}
