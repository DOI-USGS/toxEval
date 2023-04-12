#' EndPoint boxplots
#'
#' The \code{plot_tox_endpoints} function creates a set of boxplots representing EAR
#' values for each endPoint based on the selected data. A subset of data is first
#' chosen by specifying a group in the filterBy argument. The
#' filterBy argument must match one of the unique options in the category.
#' For example, if the category is "Chemical Class", then the filterBy argument
#' must be one of the defined "Chemical Class" options such as "Herbicide".

#' A boxplot is generated for each endPoint. The EAR values that are used to
#' create the boxplots are the mean or maximum (as defined by mean_logic) for each
#' site as described in "Summarizing the data"in the Introduction vignette:
#' \href{../doc/Introduction.html#summarize_data}{\code{vignette("Introduction", package = "toxEval")}}.
#'
#' Box plots are standard Tukey representations. See "Box plot details" in the Basic Workflow vignette:
#' \href{../doc/basicWorkflow.html#box_plot_details}{\code{vignette("basicWorkflow", package = "toxEval")}}
#' for more information.
#'
#' @param chemical_summary Data frame from \code{\link{get_chemical_summary}}.
#' @param category Either "Biological", "Chemical Class", or "Chemical".
#' @param filterBy Character. Either "All" or one of the filtered categories.
#' @param manual_remove Vector of categories to remove.
#' @param mean_logic Logical. \code{TRUE} displays the mean sample from each site,
#' \code{FALSE} displays the maximum sample from each site.
#' @param sum_logic logical. \code{TRUE} sums the EARs in a specified grouping,
#' \code{FALSE} does not. \code{FALSE} may be better for traditional benchmarks as
#' opposed to ToxCast benchmarks.
#' @param hit_threshold Numeric threshold defining a "hit".
#' @param font_size Numeric to adjust the axis font size.
#' @param title Character title for plot.
#' @param x_label Character for x label. Default is NA which produces an automatic label.
#' @param palette Vector of color palette for fill. Can be a named vector
#' to specify specific color for specific categories.
#' @param top_num Integer number of endpoints to include in the graph. If NA, all
#' endpoints will be included.
#' @export
#' @import ggplot2
#' @importFrom stats median
#' @examples
#' # This is the example workflow:
#' path_to_tox <- system.file("extdata", package = "toxEval")
#' file_name <- "OWC_data_fromSup.xlsx"
#'
#' full_path <- file.path(path_to_tox, file_name)
#'
#' tox_list <- create_toxEval(full_path)
#' ACC <- get_ACC(tox_list$chem_info$CAS)
#' ACC <- remove_flags(ACC)
#'
#' cleaned_ep <- clean_endPoint_info(end_point_info)
#' filtered_ep <- filter_groups(cleaned_ep)
#' chemical_summary <- get_chemical_summary(tox_list, ACC, filtered_ep)
#'
#' plot_tox_endpoints(chemical_summary,
#'   filterBy = "Cell Cycle",
#'   top_num = 10
#' )
#' \donttest{
#' plot_tox_endpoints(chemical_summary,
#'   filterBy = "Cell Cycle",
#'   top_num = 10,
#'   x_label = "EAR"
#' )
#' 
#' plot_tox_endpoints(chemical_summary,
#'   category = "Chemical Class", filterBy = "PAHs",
#'   top_num = 10, hit_threshold = 0.001
#' )
#' 
#' plot_tox_endpoints(chemical_summary, category = "Chemical", filterBy = "Atrazine")
#' plot_tox_endpoints(chemical_summary, category = "Chemical", top_num = 10)
#' single_site <- dplyr::filter(chemical_summary, site == "USGS-04024000")
#' plot_tox_endpoints(single_site, category = "Chemical", top_num = 10)
#' }
plot_tox_endpoints <- function(chemical_summary,
                               category = "Biological",
                               filterBy = "All",
                               manual_remove = NULL,
                               hit_threshold = NA,
                               mean_logic = FALSE,
                               sum_logic = TRUE,
                               font_size = NA,
                               title = NA,
                               x_label = NA,
                               palette = NA,
                               top_num = NA) {
  match.arg(category, c("Biological", "Chemical Class", "Chemical"))

  site <- endPoint <- EAR <- sumEAR <- meanEAR <- x <- y <- ".dplyr"
  CAS <- hit_label <- nonZero <- hits <- ymin <- ymax <- logEAR <- ".dplyr"

  if (nrow(chemical_summary) == 0) {
    stop("No rows in the chemical_summary data frame")
  }

  chemical_summary$EAR[chemical_summary$EAR == 0] <- NA

  if (category == "Biological") {
    chemical_summary$category <- chemical_summary$Bio_category
  } else if (category == "Chemical Class") {
    chemical_summary$category <- chemical_summary$Class
  } else {
    chemical_summary$category <- chemical_summary$chnm
  }

  single_site <- length(unique(chemical_summary$site)) == 1

  if (filterBy != "All") {
    if (!(filterBy %in% unique(chemical_summary$category))) {
      stop("filterBy argument doesn't match data")
    }
    chemical_summary <- chemical_summary[chemical_summary["category"] == filterBy, ]
  }

  if (is.na(x_label)) {
    y_label <- fancyLabels(category, mean_logic, sum_logic, single_site)
  } else {
    y_label <- x_label
  }

  if (single_site) {
    orderEP_df <- orderEP(rename(chemical_summary, meanEAR = EAR))
    orderedLevelsEP <- orderEP_df$endPoint

    if (!is.na(top_num) && top_num < length(orderedLevelsEP)) {
      orderedLevelsEP <- orderedLevelsEP[(length(orderedLevelsEP) - top_num + 1):length(orderedLevelsEP)]
      chemical_summary <- chemical_summary[chemical_summary[["endPoint"]] %in% orderedLevelsEP, ]
    }

    pretty_logs_new <- prettyLogs(chemical_summary$EAR[!is.na(chemical_summary$EAR)])

    chemical_summary$endPoint <- factor(chemical_summary$endPoint,
      levels = orderedLevelsEP
    )

    countNonZero <- chemical_summary %>%
      mutate(
        ymin = min(EAR[!is.na(EAR)], na.rm = TRUE),
        ymax = max(EAR[!is.na(EAR)], na.rm = TRUE)
      ) %>%
      group_by(endPoint, ymin, ymax) %>%
      summarise(
        nonZero = as.character(length(unique(CAS[!is.na(EAR)]))),
        hits = as.character(sum(EAR > hit_threshold, na.rm = TRUE))
      ) %>%
      ungroup()

    countNonZero$hits[countNonZero$hits == "0"] <- ""

    stackedPlot <- ggplot(data = chemical_summary) +
      theme_bw() +
      theme(
        axis.text = element_text(color = "black"),
        axis.title.y = element_blank(),
        panel.background = element_blank(),
        plot.background = element_rect(fill = "transparent", colour = NA),
        strip.background = element_rect(fill = "transparent", colour = NA),
        strip.text.y = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(vjust = .25, hjust = 1)
      ) +
      geom_hline(
        yintercept = hit_threshold, linetype = "dashed",
        color = "black", na.rm = TRUE
      )

    if (!all(is.na(palette))) {
      stackedPlot <- stackedPlot +
        geom_boxplot(aes(x = endPoint, y = EAR, fill = endPoint),
          na.rm = TRUE
        ) +
        scale_fill_manual(values = palette) +
        theme(legend.position = "none")
    } else {
      stackedPlot <- stackedPlot +
        geom_boxplot(aes(x = endPoint, y = EAR),
          fill = "steelblue", na.rm = TRUE
        )
    }
  } else {
    if (!sum_logic) {
      graphData <- chemical_summary %>%
        group_by(site, category, endPoint) %>%
        summarise(meanEAR = ifelse(mean_logic, mean(EAR, na.rm = TRUE), max(EAR, na.rm = TRUE))) %>%
        ungroup() %>%
        mutate(category = as.character(category))
    } else {
      graphData <- chemical_summary %>%
        group_by(site, date, category, endPoint) %>%
        summarise(sumEAR = sum(EAR, na.rm = TRUE)) %>%
        ungroup() %>%
        group_by(site, category, endPoint) %>%
        summarise(meanEAR = ifelse(mean_logic, mean(sumEAR, na.rm = TRUE), max(sumEAR, na.rm = TRUE))) %>%
        ungroup() %>%
        mutate(category = as.character(category))
    }

    orderEP_df <- orderEP(graphData)

    orderedLevelsEP <- orderEP_df$endPoint

    if (!is.na(top_num) &&
      top_num < length(orderedLevelsEP)) {
      orderedLevelsEP <- orderedLevelsEP[(length(orderedLevelsEP) - top_num + 1):length(orderedLevelsEP)]
      graphData <- graphData[graphData[["endPoint"]] %in% orderedLevelsEP, ]
    }

    graphData$endPoint <- factor(graphData$endPoint, levels = orderEP_df$endPoint)

    pretty_logs_new <- prettyLogs(graphData$meanEAR[!is.na(graphData$meanEAR)])
    graphData$meanEAR[graphData$meanEAR == 0] <- NA

    countNonZero <- graphData %>%
      mutate(
        ymin = min(meanEAR[!is.na(meanEAR)], na.rm = TRUE),
        ymax = max(meanEAR[!is.na(meanEAR)], na.rm = TRUE)
      ) %>%
      group_by(endPoint, ymin, ymax) %>%
      summarise(
        nonZero = as.character(length(unique(site[!is.na(meanEAR)]))),
        hits = as.character(sum(meanEAR > hit_threshold, na.rm = TRUE))
      ) %>%
      ungroup()

    countNonZero$hits[countNonZero$hits == "0"] <- ""

    stackedPlot <- ggplot(graphData) +
      theme_bw() +
      theme(
        axis.text = element_text(color = "black"),
        axis.title.y = element_blank(),
        panel.background = element_blank(),
        plot.background = element_rect(fill = "transparent", colour = NA),
        strip.background = element_rect(fill = "transparent", colour = NA),
        strip.text.y = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(vjust = .25, hjust = 1)
      )

    if (!is.na(hit_threshold)) {
      stackedPlot <- stackedPlot +
        geom_hline(
          yintercept = hit_threshold, na.rm = TRUE,
          linetype = "dashed", color = "black"
        )
    }

    if (!all(is.na(palette))) {
      stackedPlot <- stackedPlot +
        geom_boxplot(aes(x = endPoint, y = meanEAR, fill = endPoint),
          na.rm = TRUE
        ) +
        scale_fill_manual(values = palette) +
        theme(legend.position = "none")
    } else {
      stackedPlot <- stackedPlot +
        geom_boxplot(aes(x = endPoint, y = meanEAR),
          na.rm = TRUE,
          fill = "steelblue"
        )
    }
  }
  if (isTRUE(y_label == "")) {
    stackedPlot <- stackedPlot +
      scale_y_log10(
        labels = fancyNumbers,
        breaks = pretty_logs_new
      ) +
      theme(axis.title.x = element_blank())
  } else {
    stackedPlot <- stackedPlot +
      scale_y_log10(y_label,
        labels = fancyNumbers,
        breaks = pretty_logs_new
      )
  }

  plot_layout <- ggplot_build(stackedPlot)$layout

  label <- "# Sites"

  if (single_site) {
    label <- "# Chemicals"
  }

  labels_df <- countNonZero %>%
    select(-endPoint, -nonZero, -hits) %>%
    distinct() %>%
    mutate(
      x = Inf,
      label = label,
      hit_label = "# Hits"
    )


  stackedPlot <- stackedPlot +
    geom_text(
      data = countNonZero,
      aes(x = endPoint, y = ymin, label = nonZero),
      size = ifelse(is.na(font_size), 3, 0.30 * font_size),
      position = position_nudge(y = -0.05)
    ) +
    geom_text(
      data = labels_df,
      aes(x = x, y = ymin, label = label),
      size = ifelse(is.na(font_size), 3, 0.30 * font_size),
      position = position_nudge(y = -0.05)
    )

  if (!is.na(hit_threshold)) {
    stackedPlot <- stackedPlot +
      geom_text(
        data = countNonZero,
        aes(x = endPoint, y = ymax, label = hits),
        size = ifelse(is.na(font_size), 3, 0.30 * font_size),
        position = position_nudge(y = -0.05)
      ) +
      geom_text(
        data = labels_df,
        aes(x = x, y = ymax, label = hit_label),
        size = ifelse(is.na(font_size), 3, 0.30 * font_size),
        position = position_nudge(y = -0.05)
      )
  }

  if (!is.na(font_size)) {
    stackedPlot <- stackedPlot +
      theme(
        axis.text = element_text(size = font_size),
        axis.title = element_text(size = font_size)
      )
  }

  if (!is.na(title)) {
    stackedPlot <- stackedPlot +
      ggtitle(title)

    if (!is.na(font_size)) {
      stackedPlot <- stackedPlot +
        theme(plot.title = element_text(size = font_size))
    }
  }

  stackedPlot <- stackedPlot +
    coord_flip(clip = "off")

  return(stackedPlot)
}
