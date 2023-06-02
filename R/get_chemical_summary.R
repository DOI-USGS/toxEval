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

  # Getting rid of NSE warnings:
  chnm <- endPoint <- ACC_value <- Value <- `Sample Date` <- SiteID <- ".dplyr"
  EAR <- `Short Name` <- CAS <- Class <- site <- casrn <- groupCol <- ".dplyr"
  Chemical <- ".dplyr"
  
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

  chemical_summary <- full_join(ACC,
    dplyr::select(chem_data,
                   CAS, SiteID, Value, `Sample Date`),
            by = "CAS") %>%
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
      select(CAS, chnm, endPoint, site, date, EAR) %>%
      filter(endPoint %in% filtered_ep$endPoint) %>%
      left_join(select(filtered_ep, endPoint, groupCol), by = "endPoint")
  } else {
    chemical_summary <- chemical_summary %>%
      select(CAS, chnm, endPoint, site, date, EAR, groupCol)
  }

  chemical_summary <- chemical_summary %>%
    left_join(distinct(select(chem_site, site = SiteID, `Short Name`)),
      by = "site"
    ) %>%
    left_join(select(chem_info, CAS, Class), by = "CAS") %>%
    rename(
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
  chnm <- Class <- logEAR <- meanEAR <- median <- max_med <- ".dplyr"

  graphData$meanEAR[graphData$meanEAR == 0] <- NA

  orderClass_df <- graphData %>%
    mutate(logEAR = log(meanEAR)) %>%
    group_by(chnm, Class) %>%
    summarise(median = quantile(logEAR[logEAR != 0], 0.5, na.rm = TRUE)) %>%
    group_by(Class) %>%
    summarise(max_med = max(median, na.rm = TRUE)) %>%
    arrange(desc(max_med))

  return(orderClass_df)
}


orderChem <- function(graphData, orderClass_df) {
  chnm <- Class <- logEAR <- meanEAR <- median <- ".dplyr"

  graphData$meanEAR[graphData$meanEAR == 0] <- NA

  orderChem_df <- graphData %>%
    mutate(logEAR = log(meanEAR)) %>%
    group_by(chnm, Class) %>%
    summarise(median = quantile(logEAR[logEAR != 0], 0.5, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(Class = factor(Class, levels = rev(as.character(orderClass_df$Class))))

  orderChem_df$median[is.na(orderChem_df$median)] <- min(orderChem_df$median, na.rm = TRUE) - 1

  orderChem_df <- arrange(orderChem_df, Class, median)

  return(orderChem_df)
}


orderEP <- function(graphData) {
  endPoint <- logEAR <- meanEAR <- median <- ".dplyr"

  graphData$meanEAR[graphData$meanEAR == 0] <- NA

  orderEP_df <- graphData %>%
    mutate(logEAR = log(meanEAR)) %>%
    group_by(endPoint) %>%
    summarise(median = quantile(logEAR[logEAR != 0], 0.5, na.rm = TRUE)) %>%
    ungroup()

  orderEP_df$median[is.na(orderEP_df$median)] <- min(orderEP_df$median, na.rm = TRUE) - 1

  orderEP_df <- arrange(orderEP_df, median)

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
#' \strong{Flag} \tab \strong{flagsShort}\cr
#' Borderline active* \tab Borderline* \cr
#' Only highest conc above baseline, active* \tab OnlyHighest* \cr
#' Only one conc above baseline, active \tab OneAbove \cr
#' Noisy data \tab Noisy \cr
#' Hit-call potentially confounded by overfitting \tab HitCall \cr
#' Gain AC50 < lowest conc & loss AC50 < mean conc* \tab GainAC50* \cr
#' Biochemical assay with < 50\% efficacy* \tab Biochemical* \cr
#' Less than 50\% efficacy \tab LessThan50 \cr
#' AC50 less than lowest concentration tested* \tab ACCLessThan* \cr
#' GNLSmodel \tab GNLSmodel \cr
#' }
#' Asterisks indicate flags removed in the function as default.
#'
#'
#' @param ACC data frame with columns: casn, chnm, endPoint, and ACC_value
#' @param flagsShort vector of flags to to trigger REMOVAL of chemical:endPoint
#' combination. Possible values are "Borderline", "OnlyHighest", "OneAbove",
#' "Noisy", "HitCall", "GainAC50", "Biochemical","LessThan50","ACCLessThan","GNLSmodel".
#' @export
#' @examples
#' CAS <- c("121-00-6", "136-85-6", "80-05-7", "84-65-1", "5436-43-1", "126-73-8")
#' ACC <- get_ACC(CAS)
#' nrow(ACC)
#' ACC <- remove_flags(ACC)
#' nrow(ACC)
remove_flags <- function(ACC, flagsShort = c(
                           "Borderline",
                           "OnlyHighest",
                           "GainAC50",
                           "Biochemical",
                           "ACCLessThan"
                         )) {
  match.arg(flagsShort,
    c(
      "Borderline",
      "OnlyHighest",
      "OneAbove",
      "Noisy",
      "HitCall",
      "GainAC50",
      "Biochemical",
      "LessThan50",
      "ACCLessThan",
      "GNLSmodel"
    ),
    several.ok = TRUE
  )

  flags <- ".dplyr"

  flag_hits <- select(ACC, flags) %>%
    mutate(
      Borderline = grepl("Borderline active", flags),
      Noisy = grepl("Noisy data", flags),
      OneAbove = grepl("Only one conc above baseline", flags),
      OnlyHighest = grepl("Only highest conc above baseline", flags),
      Biochemical = grepl("Biochemical assay with", flags),
      GainAC50 = grepl("Gain AC50", flags),
      HitCall = grepl("potentially confounded by overfitting", flags),
      LessThan50 = grepl("Less than 50% efficacy", flags),
      ACCLessThan = grepl("AC50 less than lowest concentration tested", flags),
      GNLSmodel = grepl("Cell viability assay fit with gnls winning model", flags)
    ) %>%
    select(-flags)

  ACC <- ACC[rowSums(flag_hits[flagsShort]) == 0, ]

  return(ACC)

  # So, with the defaults, we are taking out:
  # c("Borderline active",
  #   "Only highest conc above baseline, active",
  #   "Gain AC50 < lowest conc & loss AC50 < mean conc",
  #   "Biochemical assay with < 50% efficacy")
  # We are leaving in with the defaults:
  # c("Hit-call potentially confounded by overfitting",
  #   "Only one conc above baseline, active",
  #   "AC50 less than lowest concentration tested",
  #   "Less than 50% efficacy",
  #   "Noisy data","Cell viability assay fit with gnls winning model")
}


exclude_points <- function(chemical_summary, exclusion) {
  CAS <- endPoint <- casrn <- ".dplyr"

  exclusion$CAS[exclusion$CAS == ""] <- NA
  exclusion$endPoint[exclusion$endPoint == ""] <- NA

  exclude_chem <- exclusion$CAS[is.na(exclusion$endPoint)]
  exclude_ep <- exclusion$endPoint[is.na(exclusion$CAS)]

  exclude_combo <- exclusion %>%
    filter(!is.na(CAS)) %>%
    filter(!is.na(endPoint))

  chem_filtered <- chemical_summary %>%
    filter(!(CAS %in% exclude_chem)) %>%
    filter(!(endPoint %in% exclude_ep))

  if (nrow(exclude_combo) > 0) {
    chem_filtered <- chem_filtered %>%
      anti_join(exclude_combo, by = c("CAS", "endPoint"))
  }

  return(chem_filtered)
}
