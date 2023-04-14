#' EndPoint boxplots with faceting option
#'
#' The \code{plot_tox_endpoints2} function creates a set of boxplots representing EAR
#' values for each endPoint based on the selected data. A subset of data is first
#' chosen by specifying a group in the \code{filterBy} argument. The
#' \code{filterBy} argument must match one of the unique options in the category.
#' For example, if the category is "Chemical Class", then the \code{filterBy} argument
#' must be one of the defined "Chemical Class" options such as "Herbicide".
#'
#' The difference between this function and the
#' \code{\link[toxEval]{plot_tox_endpoints}} is that
#' the \dots arguments allow for customized faceting. To include this in
#' the original toxEval function, backward compatibility would be broken.
#'
#' @param cs Data.frame from \code{\link[toxEval]{get_chemical_summary}}.
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
#' @param ... Additional group_by arguments. This can be handy for creating facet graphs.
#' @export
#' @import ggplot2
#' @importFrom stats median
#' @import dplyr
#' @examples
#' 
#' \donttest{
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
#' cs <- get_chemical_summary(tox_list, ACC, filtered_ep)
#' cs$guide_side <- "A"
#'
#' cs2 <- cs
#' cs2$guide_side <- "B"
#'
#' cs_double <- rbind(cs, cs2)
#'
#' plot_tox_endpoints2(cs_double, guide_side,
#'   top_num = 10
#' ) +
#'   ggplot2::facet_grid(. ~ guide_side, scales = "free_x")
#' }
plot_tox_endpoints2 <- function(cs, ...,
                                category = "Chemical",
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

  if (nrow(cs) == 0) {
    stop("No rows in the chemical_summary data frame")
  }

  cs$EAR[cs$EAR == 0] <- NA

  if (category == "Biological") {
    cs$category <- cs$Bio_category
  } else if (category == "Chemical Class") {
    cs$category <- cs$Class
  } else {
    cs$category <- cs$chnm
  }

  if (filterBy != "All") {
    if (!(filterBy %in% unique(cs$category))) {
      stop("filterBy argument doesn't match data")
    }
    cs <- cs[cs["category"] == filterBy, ]
  }

  if (is.na(x_label)) {
    y_label <- fancyLabels(category, mean_logic, sum_logic, FALSE)
  } else {
    y_label <- x_label
  }


  if (!sum_logic) {
    graphData <- cs %>%
      group_by(site, category, endPoint, ...) %>%
      summarise(meanEAR = ifelse(mean_logic, mean(EAR, na.rm = TRUE), max(EAR, na.rm = TRUE))) %>%
      ungroup() %>%
      mutate(category = as.character(category))
  } else {
    graphData <- cs %>%
      group_by(site, date, category, endPoint, ...) %>%
      summarise(sumEAR = sum(EAR, na.rm = TRUE)) %>%
      ungroup() %>%
      group_by(site, category, endPoint, ...) %>%
      summarise(meanEAR = ifelse(mean_logic, mean(sumEAR, na.rm = TRUE), max(sumEAR, na.rm = TRUE))) %>%
      ungroup() %>%
      mutate(category = as.character(category))
  }

  orderEP_df <- orderEP(graphData)

  orderedLevelsEP <- orderEP_df$endPoint

  if (!is.na(top_num)) {
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
    group_by(endPoint, ymin, ymax, ...) %>%
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
