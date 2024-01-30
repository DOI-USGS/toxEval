#' Prepare boxplot data
#'
#' A set of functions to prepare the data for boxplots. Often, these
#' functions are used within the plotting functions. They are exported however
#' to allow custom graphs to be created.
#'
#' The function side_by_side_data will combine two data frames,
#' either the output of \code{\link{get_chemical_summary}} or \code{\link{graph_chem_data}},
#' into a single data frame. The important work here is that the chemicals
#' and classes factor levels are ordered primarily based on "gd_left", but
#' include "gd_right" when the contents are mismatched.
#'
#' @export
#' @rdname graph_data_prep
#'
#' @param gd_left Data frame that must include the columns chnm, Class, and either EAR or meanEAR.
#' @param gd_right Data frame that must include the columns chnm, Class, and either EAR or meanEAR.
#' @param left_title Character that will be associated with the "gd_left" data
#' frame in a column named "guide_side".
#' @param right_title Character that will be associated with the "gd_right" data
#' frame in a column named "guide_side".
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
#' # Let's say we want to compare 2 chemical summaries
#' # We'll look at one summing EARs, and with concentrations
#' # First, we need a chemical summary for concentrations:
#' chemical_summary_conc <- get_concentration_summary(tox_list)
#'
#' gd_tox <- graph_chem_data(chemical_summary)
#' gd_conc <- graph_chem_data(chemical_summary_conc)
#'
#' ch_combo <- side_by_side_data(gd_tox, gd_conc,
#'   left_title = "ToxCast",
#'   right_title = "Concentrations"
#' )
#' plot_chemical_boxplots(ch_combo, guide_side,
#'   x_label = ""
#' ) +
#'   ggplot2::facet_grid(. ~ guide_side, scales = "free_x")
side_by_side_data <- function(gd_left,
                              gd_right,
                              left_title = "Left",
                              right_title = "Right") {

  gd_left$guide_side <- left_title
  gd_right$guide_side <- right_title

  gd_left_no_factor <- gd_left %>%
    dplyr::mutate(
      chnm = as.character(chnm),
      Class = as.character(Class)
    )

  gd_right_no_factor <- gd_right %>%
    dplyr::mutate(
      chnm = as.character(chnm),
      Class = as.character(Class)
    )

  chem_data_no_factors <- gd_left_no_factor %>%
    dplyr::bind_rows(gd_right_no_factor %>%
                       dplyr::filter(!(chnm %in% unique(gd_left_no_factor$chnm))))

  graph_data <- "meanEAR" %in% names(chem_data_no_factors)

  if (!(graph_data)) {
    chem_data_no_factors <- chem_data_no_factors %>%
      dplyr::rename(meanEAR = EAR)
  }

  orderChem_1_2 <- chem_data_no_factors %>%
    dplyr::group_by(chnm, Class) %>%
    dplyr::summarise(median = stats::quantile(meanEAR[meanEAR != 0], 0.5)) %>%
    dplyr::ungroup()

  class_order <- orderClass(chem_data_no_factors)

  orderChem_1_2 <- orderChem_1_2 %>%
    dplyr::mutate(Class = factor(Class, levels = class_order$Class)) %>%
    dplyr::arrange(dplyr::desc(Class), dplyr::desc(!is.na(median)), median)

  gd_1_2 <- dplyr::bind_rows(
    gd_left_no_factor,
    gd_right_no_factor
  )
  gd_1_2$Class <- factor(gd_1_2$Class, levels = class_order$Class)
  gd_1_2$chnm <- factor(gd_1_2$chnm, levels = orderChem_1_2$chnm)

  gd_1_2$guide_side <- factor(gd_1_2$guide_side,
    levels = c(left_title, right_title)
  )

  return(gd_1_2)
}

#' Create concentration summary
#'
#' Use this function to create a chemical_summary, but instead
#' of using any benchmarks, the EAR column is simply
#' the concentration. The output of this function can be used
#' in any of the plotting or table functions in the same way
#' that the output of \code{\link{get_chemical_summary}}.
#'
#'
#' @param tox_list List with data frames for chem_data, chem_info, and chem_site.
#' Created with \code{\link{create_toxEval}}.
#' @param chem_data \emph{Optional} data frame with (at least) columns: CAS, SiteID, and Value. Default is \code{NULL}.
#' The argument will over-ride what is in tox_list.
#' @param chem_site \emph{Optional} data frame with (at least) columns: SiteID, and Short Name. Default is \code{NULL}.
#' The argument will over-ride what is in tox_list.
#' @param chem_info \emph{Optional} data frame with (at least) columns: CAS, and class. Default is \code{NULL}.
#' The argument will over-ride what is in tox_list.
#' @param tox_names Logical whether to use the provided chemical names from the ToxCast or not. If
#' there is not a match by CAS, the function will look for a column "Chemical" in the "Chemical"
#' tab. If that column doesn't exist, it will create a (not good!) name.
#' @export
#' @return a data frame with the columns: CAS, chnm (chemical name
#' as a factor), site, date, EAR (which is just concentration), Bio_category, shortName (of site), Class. The output of this
#' function is where you find EAR values for every chemical/endpoint combination.
#'
#' @examples
#' path_to_tox <- system.file("extdata", package = "toxEval")
#' file_name <- "OWC_data_fromSup.xlsx"
#' full_path <- file.path(path_to_tox, file_name)
#'
#' tox_list <- create_toxEval(full_path)
#'
#' chemical_summary_conc <- get_concentration_summary(tox_list)
#' head(chemical_summary_conc)
#' plot_tox_boxplots(chemical_summary_conc,
#'   category = "Chemical",
#'   x_label = "Concentration [ug/L]"
#' )
get_concentration_summary <- function(tox_list,
                                      chem_data = NULL,
                                      chem_site = NULL,
                                      chem_info = NULL,
                                      tox_names = FALSE) {

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

  if (is.character(chem_data$Value)) {
    chem_data$Value <- as.numeric(chem_data$Value)
  }

  chemical_summary <- chem_data %>%
    dplyr::select(CAS, SiteID, Value, `Sample Date`) %>%
    dplyr::filter(!is.na(Value)) %>%
    dplyr::rename(
      EAR = Value,
      site = SiteID,
      date = `Sample Date`
    ) %>%
    dplyr::mutate(
      Bio_category = "Concentration",
      endPoint = "Concentration"
    ) %>%
    dplyr::left_join(dplyr::select(chem_site,
                     site = SiteID,
                     shortName = `Short Name`),
      by = "site") %>%
    dplyr::left_join(dplyr::select(chem_info, CAS, Class, chnm = Chemical),
              by = "CAS")

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
