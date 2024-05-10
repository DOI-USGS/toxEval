#' Get the ACC values for a selection of chemicals
#'
#' The \code{get_ACC} function retrieves the activity concentration at cutoff
#' (ACC) values for specified chemicals.
#'
#' The function \code{get_ACC} will
#' convert the ACC values in the ToxCast database from units of (\eqn{\mu}M)
#' to units of \eqn{\mu}g/L, and reformat the data as input to toxEval.
#'
#' @param CAS Vector of CAS.
#'
#' @return data frame with columns CAS, chnm, flags, endPoint, ACC, MlWt, and ACC_value
#' @export
#' @examples
#' CAS <- c("121-00-6", "136-85-6", "80-05-7", "84-65-1", "5436-43-1", "126-73-8")
#' ACC <- get_ACC(CAS)
#' head(ACC)
get_ACC <- function(CAS) {

  chem_list <- tox_chemicals |> 
    dplyr::select(casrn = .data$casn,
                  MlWt = .data$Structure_MolWt) |> 
    dplyr::filter(.data$casrn %in% CAS)

  ACC <- ToxCast_ACC |> 
    dplyr::filter(.data$casn %in% CAS) |> 
    dplyr::right_join(chem_list, by = c("casn" = "casrn")) |> 
    dplyr::rename(CAS = casn) |> 
    dplyr::mutate(ACC_value = .data$hit_val * .data$MlWt) |> 
    dplyr::filter(!is.na(.data$ACC_value)) |> 
    dplyr::left_join(dplyr::select(tox_chemicals,
                                  CAS = .data$casn,
                                  chnm = .data$chnm),
                     by = "CAS") |> 
    dplyr::left_join(end_point_info |> 
                       dplyr::select(.data$aeid,
                                     endPoint = .data$assay_component_endpoint_name),
                     by = "aeid")

  if (any(is.na(ACC$MlWt))) {
    warning("Some chemicals are missing molecular weights")
  }

  return(ACC)
}
