#' Rank endpoints by category
#'
#' The \code{endpoint_hits_DT} (data.table (DT) option) and \code{endpoint_hits}
#' (data frame option) functions create tables with one row per endPoint, and
#' one column per category("Biological", "Chemical", or "Chemical Class"). The
#' values in the table are the number of sites where the EAR exceeded the
#' user-defined EAR hit_threshold in that endpoint/category combination. If the
#' category "Chemical" is chosen, an "info" link is provided to the
#' chemical information available in the "Comptox Dashboard"
#' \url{https://comptox.epa.gov/dashboard/}.
#'
#' The tables show slightly different results when choosing to explore data
#' from a single site rather than all sites. The value displayed in this
#' instance is the number of samples with hits rather than the number of sites
#' with hits.
#'
#' @param chemical_summary Data frame from \code{get_chemical_summary}
#' @param mean_logic Logical.  \code{TRUE} displays the mean sample from each site,
#' FALSE displays the maximum sample from each site.
#' @param sum_logic Logical. \code{TRUE} sums the EARs in a specified grouping,
#' \code{FALSE} does not. \code{FALSE} indicates that EAR values are not considered to be
#' additive and often will be a more appropriate choice for traditional
#' benchmarks as opposed to ToxCast benchmarks.
#' @param category Character. Either "Biological", "Chemical Class", or "Chemical".
#' @param hit_threshold Numeric. EAR threshold defining a "hit".
#' @param include_links Logical. whether or not to include a link to the ToxCast
#' dashboard. Only needed for the "Chemical" category.
#' @export
#' @return data frame with one row per endpoint that had a hit (based on the
#' hit_threshold). The columns are based on the category.
#' @rdname endpoint_hits_DT
#' @importFrom stats median
#' @examples
#' # This is the example workflow:
#' path_to_tox <- system.file("extdata", package = "toxEval")
#' file_name <- "OWC_data_fromSup.xlsx"
#'
#' full_path <- file.path(path_to_tox, file_name)
#' \donttest{
#' tox_list <- create_toxEval(full_path)
#'
#' ACC <- get_ACC(tox_list$chem_info$CAS)
#' ACC <- remove_flags(ACC)
#'
#' cleaned_ep <- clean_endPoint_info(end_point_info)
#' filtered_ep <- filter_groups(cleaned_ep)
#' chemical_summary <- get_chemical_summary(tox_list, ACC, filtered_ep)
#'
#' hits_df <- endpoint_hits(chemical_summary, category = "Biological")
#' endpoint_hits_DT(chemical_summary, category = "Biological")
#' endpoint_hits_DT(chemical_summary, category = "Chemical Class")
#' endpoint_hits_DT(chemical_summary, category = "Chemical")
#' }
endpoint_hits_DT <- function(chemical_summary,
                             category = "Biological",
                             mean_logic = FALSE,
                             sum_logic = TRUE,
                             hit_threshold = 0.1,
                             include_links = TRUE) {

  fullData <- endpoint_hits(
    chemical_summary = chemical_summary,
    category = category,
    mean_logic = mean_logic,
    sum_logic = sum_logic,
    hit_threshold = hit_threshold
  )

  if (category == "Chemical") {
    orig_names <- names(fullData)

    casKey <- dplyr::select(chemical_summary, chnm, CAS) %>%
      dplyr::distinct()

    numeric_hits <- fullData
    hits <- sapply(fullData, function(x) as.character(x))

    if (include_links) {
      for (k in 1:nrow(fullData)) {
        for (z in 2:ncol(fullData)) {
          if (!is.na(fullData[k, z])) {
            if (fullData[k, z] < 10) {
              hit_char <- paste0("0", fullData[k, z])
            } else {
              hit_char <- as.character(fullData[k, z])
            }
            hits[k, z] <- paste(hit_char, createLink(cas = casKey$CAS[casKey$chnm == names(fullData)[z]]))
          }
        }
      }
    }
    fullData <- data.frame(hits, stringsAsFactors = FALSE)
    names(fullData) <- orig_names
  }

  n <- ncol(fullData) - 1

  if (n > 20 & n < 30) {
    colors <- c(
      RColorBrewer::brewer.pal(n = 12, name = "Set3"),
      RColorBrewer::brewer.pal(n = 8, name = "Set2"),
      RColorBrewer::brewer.pal(n = max(c(3, n - 20)), name = "Set1")
    )
  } else if (n <= 20) {
    colors <- c(
      RColorBrewer::brewer.pal(n = 12, name = "Set3"),
      RColorBrewer::brewer.pal(n = max(c(3, n - 12)), name = "Set2")
    )
  } else {
    colors <- colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))(n)
  }

  fullData_dt <- DT::datatable(fullData,
    extensions = "Buttons",
    escape = FALSE,
    rownames = FALSE,
    options = list(
      dom = "Bfrtip",
      buttons = list("colvis"),
      scrollX = TRUE,
      order = list(list(1, "desc"))
    )
  )

  for (i in 2:ncol(fullData)) {
    fullData_dt <- DT::formatStyle(fullData_dt,
      names(fullData)[i],
      backgroundColor = colors[i]
    )

    if (category != "Chemical") {
      fullData_dt <- DT::formatStyle(fullData_dt, names(fullData)[i],
        background = DT::styleColorBar(range(fullData[, names(fullData)[i]], na.rm = TRUE), "goldenrod"),
        backgroundSize = "100% 90%",
        backgroundRepeat = "no-repeat",
        backgroundPosition = "center"
      )
    }
  }

  return(fullData_dt)
}

#' @rdname endpoint_hits_DT
#' @export
endpoint_hits <- function(chemical_summary,
                          category = "Biological",
                          mean_logic = FALSE,
                          sum_logic = TRUE,
                          hit_threshold = 0.1) {

  match.arg(category, c("Biological", "Chemical Class", "Chemical"))

  fullData_init <- data.frame(endPoint = "", stringsAsFactors = FALSE)
  fullData <- fullData_init

  if (category == "Chemical") {
    chemical_summary <- dplyr::mutate(chemical_summary, category = chnm)
  } else if (category == "Chemical Class") {
    chemical_summary <- dplyr::mutate(chemical_summary, category = Class)
  } else {
    chemical_summary <- dplyr::mutate(chemical_summary, category = Bio_category)
  }

  if (length(unique(chemical_summary$site)) > 1) {
    if (!sum_logic) {
      fullData <- chemical_summary %>%
        dplyr::group_by(site, category, endPoint, date) %>%
        dplyr::summarize(sumEAR = max(EAR)) %>%
        dplyr::group_by(site, category, endPoint) %>%
        dplyr::summarize(meanEAR = ifelse(mean_logic, mean(sumEAR), max(sumEAR))) %>%
        dplyr::group_by(category, endPoint) %>%
        dplyr::summarize(nSites = sum(meanEAR > hit_threshold)) %>%
        tidyr::spread(category, nSites)
    } else {
      fullData <- chemical_summary %>%
        dplyr::group_by(site, category, endPoint, date) %>%
        dplyr::summarize(sumEAR = sum(EAR)) %>%
        dplyr::group_by(site, category, endPoint) %>%
        dplyr::summarize(meanEAR = ifelse(mean_logic, mean(sumEAR), max(sumEAR))) %>%
        dplyr::group_by(category, endPoint) %>%
        dplyr::summarize(nSites = sum(meanEAR > hit_threshold)) %>%
        tidyr::spread(category, nSites)
    }
  } else {
    if (!sum_logic) {
      fullData <- chemical_summary %>%
        dplyr::group_by(category, endPoint) %>%
        dplyr::summarise(nSites = sum(EAR > hit_threshold)) %>%
        tidyr::spread(category, nSites)
    } else {
      fullData <- chemical_summary %>%
        dplyr::group_by(category, endPoint, date) %>%
        dplyr::summarize(sumEAR = sum(EAR)) %>%
        dplyr::group_by(category, endPoint) %>%
        dplyr::summarise(nSites = sum(sumEAR > hit_threshold)) %>%
        tidyr::spread(category, nSites)
    }
  }

  if (any(rowSums(fullData[, -1], na.rm = TRUE) > 0)) {
    fullData <- fullData[(rowSums(fullData[, -1], na.rm = TRUE) != 0), ]
  }

  fullData <- fullData[, colSums(is.na(fullData)) != nrow(fullData)]

  sumOfColumns <- colSums(fullData[c(-1)], na.rm = TRUE)
  if (!all(sumOfColumns == 0)) {
    orderData <- order(sumOfColumns, decreasing = TRUE)
    orderData <- orderData[sumOfColumns[orderData] != 0] + 1

    fullData <- fullData[, c(1, orderData)]
  }

  fullData <- fullData[order(fullData[[2]], decreasing = TRUE), ]

  return(fullData)
}

#' createLink
#'
#' Create links
#' @param cas character
#' @export
#' @keywords internal
createLink <- function(cas) {
  paste0('<a href="https://comptox.epa.gov/dashboard/dsstoxdb/results?search=', cas, '#invitrodb-bioassays-toxcast-tox21" target="_blank">&#9432;</a>')
}
