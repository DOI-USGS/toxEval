#' Biological hits per category
#'
#' The \code{hits_by_groupings_DT} (DT option) and
#' \code{hits_by_groupings} (data frame option) functions create tables
#' with one row per category("Biological", "Chemical", or "Chemical Class").
#' The columns indicate the "Biological" groupings. The values in the table
#' signify how many sites have samples with EARs that exceeded the hit_threshold
#' for that particular "Biological"/category combination. If the user chooses
#' "Biological" as the category, it is a simple 2-column table of "Biological"
#' groupings and number of sites (nSites).
#'
#' The tables result in slightly different results for a single site, displaying
#' the number of samples with hits rather than the number of sites.
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
#' @return data frame with one row per category, and one column per Biological grouping.
#' @rdname hits_by_groupings_DT
#' @examples
#' # This is the example workflow:
#' path_to_tox <- system.file("extdata", package = "toxEval")
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
#' site_df <- hits_by_groupings(chemical_summary, category = "Biological")
#'
#' hits_by_groupings_DT(chemical_summary, category = "Biological")
#' hits_by_groupings_DT(chemical_summary, category = "Chemical Class")
#' hits_by_groupings_DT(chemical_summary, category = "Chemical")
#'
hits_by_groupings_DT <- function(chemical_summary,
                                 category = "Biological",
                                 mean_logic = FALSE,
                                 sum_logic = TRUE,
                                 hit_threshold = 0.1) {
  match.arg(category, c("Biological", "Chemical Class", "Chemical"))

  tableData <- hits_by_groupings(
    chemical_summary = chemical_summary,
    category = category,
    mean_logic = mean_logic,
    sum_logic = sum_logic,
    hit_threshold = hit_threshold
  )

  cuts <- seq(0, max(as.matrix(tableData), na.rm = TRUE), length.out = 8)
  colors <- RColorBrewer::brewer.pal(9, "Blues") # "RdYlBu"

  tableData1 <- DT::datatable(tableData,
    extensions = "Buttons",
    rownames = TRUE,
    options = list(
      scrollX = TRUE,
      dom = "Bfrtip",
      buttons = list("colvis"),
      order = list(list(1, "desc"))
    )
  )

  if (category != "Biological") {
    for (i in 1:ncol(tableData)) {
      tableData1 <- DT::formatStyle(tableData1,
        columns = names(tableData)[i],
        backgroundColor = DT::styleInterval(cuts = cuts, values = colors),
        color = DT::styleInterval(0.75 * max(tableData, na.rm = TRUE), values = c("black", "white")),
        `font-size` = "17px"
      )
    }
  }

  return(tableData1)
}

#' @export
#' @rdname hits_by_groupings_DT
hits_by_groupings <- function(chemical_summary,
                              category,
                              mean_logic = FALSE,
                              sum_logic = TRUE,
                              hit_threshold = 0.1) {

  match.arg(category, c("Biological", "Chemical Class", "Chemical"))

  if (category == "Biological") {
    chemical_summary$category <- chemical_summary$Bio_category
  } else if (category == "Chemical Class") {
    chemical_summary$category <- chemical_summary$Class
  } else {
    chemical_summary$category <- chemical_summary$chnm
  }

  if (length(unique(chemical_summary$site)) > 1) {
    if (!sum_logic) {
      tableData <- chemical_summary %>%
        dplyr::group_by(site, Bio_category, category) %>%
        dplyr::summarize(meanEAR = ifelse(mean_logic, mean(EAR), max(EAR))) %>%
        dplyr::group_by(Bio_category, category) %>%
        dplyr::summarize(nSites = sum(meanEAR > hit_threshold)) %>%
        data.frame()
    } else {
      tableData <- chemical_summary %>%
        dplyr::group_by(site, Bio_category, category, date) %>%
        dplyr::summarize(sumEAR = sum(EAR)) %>%
        dplyr::group_by(site, Bio_category, category) %>%
        dplyr::summarize(meanEAR = ifelse(mean_logic, mean(sumEAR), max(sumEAR))) %>%
        dplyr::group_by(Bio_category, category) %>%
        dplyr::summarize(nSites = sum(meanEAR > hit_threshold)) %>%
        data.frame()
    }
  } else {
    if (!sum_logic) {
      tableData <- chemical_summary %>%
        dplyr::group_by(Bio_category, category) %>%
        dplyr::summarise(nSites = sum(EAR > hit_threshold)) %>%
        data.frame()
    } else {
      tableData <- chemical_summary %>%
        dplyr::group_by(Bio_category, category, date) %>%
        dplyr::summarise(sumEAR = sum(EAR)) %>%
        data.frame() %>%
        dplyr::group_by(Bio_category, category) %>%
        dplyr::summarise(nSites = sum(sumEAR > hit_threshold)) %>%
        data.frame()
    }
  }

  if (category != "Biological") {
    tableData <- tableData %>%
      tidyr::spread(Bio_category, nSites)

    sumOfColumns <- colSums(tableData[-1], na.rm = TRUE)

    if (!all(sumOfColumns == 0)) {
      orderData <- order(sumOfColumns, decreasing = TRUE)
      orderData <- orderData[sumOfColumns[orderData] != 0] + 1

      tableData <- tableData[, c(1, orderData)]
    }

    groups <- tableData$category

    tableData <- tableData[!is.na(groups), -1, drop = FALSE]
    rownames(tableData) <- groups[!is.na(groups)]
  } else {
    tableData <- dplyr::select(tableData, Bio_category, nSites)
    rownames(tableData) <- tableData$Bio_category
    tableData <- tableData[, -1, drop = FALSE]
  }

  if (length(unique(chemical_summary$site)) == 1) {
    names(tableData)[names(tableData) == "nSites"] <- "nSamples"
  }

  return(tableData)
}
