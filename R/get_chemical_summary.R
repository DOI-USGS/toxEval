#' Compute EAR values
#'
#' This function computes Exposure:Activity ratios using user-provided measured
#' concentration data from the output of \code{\link{create_toxEval}},
#' and joins the data with the activity concentration at cutoff data provided by
#' ToxCast.Data from ToxCast is included with this package, but alternative
#' benchmark data can be provided to perform the same "toxEval" analysis.
#'
#' To use the data provided by the package, a sample workflow is shown below
#' in the examples. The examples include retrieving the ToxCast (ACC) values
#' that are used to calculate EARs, choosing endPoints that should be ignored
#' based on data quality "flags" in the ToxCast database, and removing groups of
#' endPoints that may not be important to the analysis at hand.
#'
#'
#' @param tox_list List with data frames for chem_data, chem_info, chem_site,
#' and optionally exclusions and benchmarks. Created with \code{\link{create_toxEval}}.
#' @param ACC Data frame with columns: CAS, chnm, endPoint, and ACC_value
#' for specific chemical/endpoint combinations generated using the
#' \code{\link{get_ACC}} function. EndPoints with specific data quality flags
#' may optionally be removed using the \code{\link{remove_flags}} function.
#' @param filtered_ep Data frame with columns: endPoints, groupCol. Default is \code{"All"}, where no
#' filtering occurs.
#' @param chem_data \emph{Optional} data frame with (at least) columns: CAS, SiteID, and Value. Default is \code{NULL}.
#' The argument will over-ride what is in tox_list.
#' @param chem_site \emph{Optional} data frame with (at least) columns: SiteID, and Short Name. Default is \code{NULL}.
#' The argument will over-ride what is in tox_list.
#' @param chem_info \emph{Optional} data frame with (at least) columns: CAS, and class. Default is \code{NULL}.
#' The argument will over-ride what is in tox_list.
#' @param exclusion \emph{Optional} data frame with (at least) columns: CAS and endPoint. Default is \code{NULL}.
#' The argument will over-ride what is in tox_list.
#' @export
#' @return a data frame with the columns: CAS, chnm (chemical name
#' as a factor), site, date, EAR, Bio_category, shortName (of site), Class. The output of this
#' function is where you find EAR values for every chemical/endpoint combination.
#'
#' @examples
#' path_to_tox <- system.file("extdata", package = "toxEval")
#' file_name <- "OWC_data_fromSup.xlsx"
#' full_path <- file.path(path_to_tox, file_name)
#'
#' tox_list <- create_toxEval(full_path)
#'
#' ACC <- get_ACC(tox_list$chem_info$CAS)
#' ACC <- remove_flags(ACC)
#'
#' cleaned_ep <- clean_endPoint_info(end_point_info)
#' filtered_ep <- filter_groups(cleaned_ep)
#'
#' chemical_summary <- get_chemical_summary(tox_list, ACC, filtered_ep)
#' head(chemical_summary)
get_chemical_summary <- function(tox_list, ACC = NULL, filtered_ep = "All",
                                 chem_data = NULL, chem_site = NULL,
                                 chem_info = NULL, exclusion = NULL) {

  if (is.null(chem_data)) {
    chem_data <- tox_list[["chem_data"]]
  } else {
    chem_data <- rm_em_dash(chem_data)
  }

  if (is.null(chem_site)) {
    chem_site <- tox_list[["chem_site"]]
  } else {
    chem_site <- rm_em_dash(chem_site)
  }

  if (is.null(chem_info)) {
    chem_info <- tox_list[["chem_info"]]
  } else {
    chem_info <- rm_em_dash(chem_info)
  }

  if (is.null(exclusion)) {
    exclusion <- tox_list[["exclusions"]]
  } else {
    exclusion <- rm_em_dash(exclusion)
  }

  if (is.null(ACC)) {
    ACC <- tox_list[["benchmarks"]]
    ACC <- dplyr::select(ACC, CAS, endPoint, ACC_value, groupCol)
  } else {
    ACC <- dplyr::select(ACC, CAS, endPoint, ACC_value)
  }
  
  
  
  if (is.character(chem_data$Value)) {
    chem_data$Value <- as.numeric(chem_data$Value)
  }

  chemical_summary <- dplyr::full_join(dplyr::distinct(ACC),
    dplyr::select(chem_data,
                   CAS, SiteID, Value, `Sample Date`) %>% 
      dplyr::filter(!is.na(CAS)),
            by = "CAS", 
            relationship = "many-to-many") %>%
    dplyr::filter(
      !is.na(ACC_value),
      !is.na(Value)) %>%
    dplyr::mutate(EAR = Value / ACC_value) %>%
    dplyr::rename(site = SiteID,
                  date = `Sample Date`) %>%
    dplyr::left_join(dplyr::select(chem_info, CAS, chnm = Chemical),
                     by = "CAS")

  if (all(filtered_ep != "All")) {
    chemical_summary <- chemical_summary %>%
      dplyr::select(CAS, chnm, endPoint, site, date, EAR) %>%
      dplyr::filter(endPoint %in% filtered_ep$endPoint) %>%
      dplyr::left_join(dplyr::select(filtered_ep, endPoint, groupCol), by = "endPoint")
  } else {
    chemical_summary <- chemical_summary %>%
      dplyr::select(CAS, chnm, endPoint, site, date, EAR, groupCol)
  }

  chemical_summary <- chemical_summary %>%
    dplyr::left_join(dplyr::distinct(dplyr::select(chem_site, site = SiteID, `Short Name`)),
      by = "site"
    ) %>%
    dplyr::left_join(dplyr::distinct(dplyr::select(chem_info, CAS, Class)), by = "CAS") %>%
    dplyr::rename(
      Bio_category = groupCol,
      shortName = `Short Name`
    )

  if (!is.null(exclusion)) {
    chemical_summary <- exclude_points(chemical_summary, exclusion)
  }

  graphData <- graph_chem_data(chemical_summary)

  orderClass_df <- orderClass(graphData)

  orderChem_df <- orderChem(graphData, orderClass_df)

  chemical_summary$chnm <- factor(chemical_summary$chnm,
    levels = unique(orderChem_df$chnm)
  )

  chemical_summary$Class <- factor(chemical_summary$Class,
    levels = rev(levels(orderChem_df$Class))
  )

  return(chemical_summary)
}


orderClass <- function(graphData) {

  graphData$meanEAR[graphData$meanEAR == 0] <- NA

  orderClass_df <- graphData %>%
    dplyr::mutate(logEAR = log(meanEAR)) %>%
    dplyr::group_by(chnm, Class) %>%
    dplyr::summarise(median = stats::quantile(logEAR[logEAR != 0], 0.5, na.rm = TRUE)) %>%
    dplyr::group_by(Class) %>%
    dplyr::summarise(max_med = max(median, na.rm = TRUE)) %>%
    dplyr::arrange(dplyr::desc(max_med))

  return(orderClass_df)
}


orderChem <- function(graphData, orderClass_df) {
  
  graphData$meanEAR[graphData$meanEAR == 0] <- NA

  orderChem_df <- graphData %>%
    dplyr::mutate(logEAR = log(meanEAR)) %>%
    dplyr::group_by(chnm, Class) %>%
    dplyr::summarise(median = stats::quantile(logEAR[logEAR != 0], 0.5, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(Class = factor(Class, levels = rev(as.character(orderClass_df$Class))))

  orderChem_df$median[is.na(orderChem_df$median)] <- min(orderChem_df$median, na.rm = TRUE) - 1

  orderChem_df <- dplyr::arrange(orderChem_df, Class, median)

  return(orderChem_df)
}


orderEP <- function(graphData) {
  
  graphData$meanEAR[graphData$meanEAR == 0] <- NA

  orderEP_df <- graphData %>%
    dplyr::mutate(logEAR = log(meanEAR)) %>%
    dplyr::group_by(endPoint) %>%
    dplyr::summarise(median = stats::quantile(logEAR[logEAR != 0], 0.5, na.rm = TRUE)) %>%
    dplyr::ungroup()

  orderEP_df$median[is.na(orderEP_df$median)] <- min(orderEP_df$median, na.rm = TRUE) - 1

  orderEP_df <- dplyr::arrange(orderEP_df, median)

  return(orderEP_df)
}

#' Remove endpoints with specific data quality flags from data
#'
#' Through the ToxCast program quality assurance procedures, information
#' is examined and at times, it is necessary to assign a data quality flag
#' to a specific chemical/assay result. A toxEval user may want to include
#' or exclude assay results with certain flags depending on the objectives
#' of a given study. Assay results with specific data quality flags assigned
#' to them can be removed based on their designated flag with the
#' \code{remove_flags} function. The flags included in ToxCast, and the associated
#' flagsShort value (used in the remove_flags function) are as follows:
#' \tabular{ll}{
#' \strong{flag_id} \tab \strong{Full Name}\cr
#'5* \tab Model directionality questionable \cr
#'6* \tab Only highest conc above baseline, active \cr
#'7 \tab Only one conc above baseline, active \cr
#'8 \tab Multiple points above baseline, inactive \cr
#'9 \tab Bmd > ac50, indication of high baseline variability \cr
#'10 \tab Noisy data \cr
#'11* \tab Borderline \cr
#'15* \tab Gain AC50 < lowest conc & loss AC50 < mean conc \cr
#'17 \tab Less than 50\% efficacy \cr
#'18* \tab AC50 less than lowest concentration tested \cr
#'13 \tab Average number of replicates per conc is less than 2 \cr
#'14 \tab Number of concentrations tested is less than 4 \cr
#'19 \tab Cell viability assay fit with gnls winning model \cr
#' }
#' Asterisks indicate flags removed in the function as default.
#'
#'
#' @param ACC data frame with columns: casn, chnm, endPoint, and ACC_value
#' @param flag_id vector of flags to to trigger REMOVAL 
#' @export
#' @examples
#' CAS <- c("121-00-6", "136-85-6", "80-05-7", "84-65-1", "5436-43-1", "126-73-8")
#' ACC <- get_ACC(CAS)
#' nrow(ACC)
#' 
#' # See available flags and associated ids:
#' 
#' flags
#' 
#' ACC <- remove_flags(ACC)
#' nrow(ACC)
remove_flags <- function(ACC, 
                         flag_id = c(5, 6, 11, 15, 18)) {
  match.arg(as.character(flag_id),
            as.character(unique(flags$flag_id)),
    several.ok = TRUE
  )

  remove_rows <- which(colSums(sapply(ACC$flags, "%in%", x = flag_id)) > 0)
  
  ACC_filter <- ACC[-remove_rows, ]

  return(ACC_filter)


}


exclude_points <- function(chemical_summary, exclusion) {

  exclusion$CAS[exclusion$CAS == ""] <- NA
  exclusion$endPoint[exclusion$endPoint == ""] <- NA

  exclude_chem <- exclusion$CAS[is.na(exclusion$endPoint)]
  exclude_ep <- exclusion$endPoint[is.na(exclusion$CAS)]

  exclude_combo <- exclusion %>%
    dplyr::filter(!is.na(CAS)) %>%
    dplyr::filter(!is.na(endPoint))

  chem_filtered <- chemical_summary %>%
    dplyr::filter(!(CAS %in% exclude_chem)) %>%
    dplyr::filter(!(endPoint %in% exclude_ep))

  if (nrow(exclude_combo) > 0) {
    chem_filtered <- chem_filtered %>%
      dplyr::anti_join(exclude_combo, by = c("CAS", "endPoint"))
  }

  return(chem_filtered)
}
