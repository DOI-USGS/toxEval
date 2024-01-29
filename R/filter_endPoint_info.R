#' Filter endPoints based on groups and assays.
#'
#' This function provides a mechanism to specify 3 levels of information in the
#' supplied data frame \code{\link{end_point_info}} to be used in subsequent analysis steps.
#' First, the user specifies the ToxCast assay annotation using the 'groupCol'
#' argument, which is a column header in 'end_point_info'. Second, the user
#' specifies the families of assays to use. Finally, the user can choose to
#' remove specific group(s) from the category. The default is to remove
#' 'Background Measurement' and 'Undefined'. Choices for this should be
#' reconsidered based on individual study objectives.
#'
#' The default category ('groupCol') is 'intended_target_family'. Depending
#' on the study, other categories may be more relevant. The best resource on these
#' groupings is the "ToxCast Assay Annotation Data User Guide".
#' It defines "intended_target_family" as "the target family of the
#' objective target for the assay". Much more detail can be discovered in that documentation.
#'
#' @param ep Data frame containing Endpoint information from ToxCast
#' @param groupCol Character name of ToxCast annotation column to use as a group category
#' @param assays Vector of assays to use in the data analysis. Possible values are "ACEA", "APR", "ATG",
#' "NVS", "OT", "TOX21", "CEETOX", "LTEA", "CLD", "TANGUAY", "CCTE_PADILLA", "BSK" ,
#' "CCTE", "STM", "ARUNA", "CCTE_SHAFER", "CPHEA_STOKER", "CCTE_GLTED", "UPITT", "UKN",
#' "ERF", "TAMU", "IUF", "CCTE_MUNDY", "UTOR", "VALA". By default, the
#' "BSK" (BioSeek) assay is removed.
#' @param remove_groups Vector of groups within the selected 'groupCol' to remove.
#' @export
#' @examples
#' end_point_info <- end_point_info
#' cleaned_ep <- clean_endPoint_info(end_point_info)
#' filtered_ep <- filter_groups(cleaned_ep)
#' head(filtered_ep)
filter_groups <- function(ep,
                          groupCol = "intended_target_family",
                          assays = c(
                            "ACEA", "APR", "ATG",
                            "NVS", "OT", "TOX21", "CEETOX",
                            "LTEA", "CLD", "TANGUAY", "CCTE_PADILLA",
                            "CCTE", "STM", "ARUNA", "CCTE_SHAFER",
                            "CPHEA_STOKER", "CCTE_GLTED", "UPITT", "UKN",
                            "ERF", "TAMU", "IUF", "CCTE_MUNDY",
                            "UTOR", "VALA"
                          ),
                          remove_groups = c("Background Measurement", "Undefined")) {
  possible_assays <- unique(end_point_info$assay_source_name)
  match.arg(assays, possible_assays, several.ok = TRUE)

  # Getting rid of NSE warnings:
  assay_source_name <- assay_component_endpoint_name <- ".dplyr"

  ep <- ep[, c("assay_component_endpoint_name", groupCol, "assay_source_name")] %>%
    dplyr::rename(
      endPoint = assay_component_endpoint_name,
      assaysFull = assay_source_name
    )
  names(ep)[names(ep) == groupCol] <- "groupCol"

  ep <- ep[(ep$assaysFull %in% assays), ]
  ep <- ep[!is.na(ep$groupCol), ]
  if (any(!is.na(remove_groups))) {
    if (any(remove_groups != "")) {
      ep <- ep[!(ep$groupCol) %in% remove_groups, ]
    }
  }

  return(ep)
}
