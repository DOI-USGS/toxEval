#' clean_endPoint_info
#'
#' Define a subset of the ToxCast database for relevance to toxEval analyses.
#' Subsetting is done based upon methods defined by
#' \doi{10.1021/acs.est.7b01613}{Blackwell et al., 2017}.
#' Specifically, this function removes endPoints that are ATG sources with
#' signal loss, and NVS with signal gain (basically: some assay/signal combinations
#' are removed because they target non-specific endpoints). Also, this function adds additional
#' categories to intended_target_family and intended_target_family_sub as
#' described in the paper linked above.
#'
#' @param end_point_info Data frame Endpoint information from ToxCast.
#' @export
#' @return The returned data frame is based on end_point_info, but with some endPoints
#' filtered out and some additional categories in intended_target_family and
#' intended_target_family_sub. The names in intended_target_family
#' are revised to look more appealing in graphs and tables.
#' @examples
#' end_point_info <- end_point_info
#' nrow(end_point_info)
#' cleaned_ep <- clean_endPoint_info(end_point_info)
#' nrow(cleaned_ep)
clean_endPoint_info <- function(end_point_info) {
  end_point_info <- end_point_info[!(end_point_info$assay_source_name == "ATG" & end_point_info$signal_direction == "loss"), ]
  end_point_info <- end_point_info[!(end_point_info$assay_source_name == "NVS" & end_point_info$signal_direction == "gain"), ]

  end_point_info$intended_target_family[end_point_info$assay_component_endpoint_name %in%
    c(
      "CLD_CYP1A1_24hr", "CLD_CYP1A1_48hr", "CLD_CYP1A1_6hr",
      "CLD_CYP1A2_24hr", "CLD_CYP1A2_48hr", "CLD_CYP1A2_6hr"
    )] <- "dna binding"

  end_point_info$intended_target_family[end_point_info$assay_component_endpoint_name %in%
    c(
      "CLD_CYP2B6_24hr", "CLD_CYP2B6_48hr", "CLD_CYP2B6_6hr",
      "CLD_CYP3A4_24hr", "CLD_CYP3A4_48hr", "CLD_CYP3A4_6hr",
      "CLD_SULT2A_48hr", "CLD_UGT1A1_48hr", "NVS_NR_bER",
      "NVS_NR_bPR", "NVS_NR_cAR"
    )] <- "nuclear receptor"

  end_point_info$intended_target_family[end_point_info$assay_component_endpoint_name %in%
    c(
      "Tanguay_ZF_120hpf_ActivityScore",
      "Tanguay_ZF_120hpf_AXIS_up",
      "Tanguay_ZF_120hpf_BRAI_up",
      "Tanguay_ZF_120hpf_CFIN_up",
      "Tanguay_ZF_120hpf_EYE_up",
      "Tanguay_ZF_120hpf_JAW_up",
      "Tanguay_ZF_120hpf_MORT_up",
      "Tanguay_ZF_120hpf_OTIC_up",
      "Tanguay_ZF_120hpf_PE_up",
      "Tanguay_ZF_120hpf_PFIN_up",
      "Tanguay_ZF_120hpf_PIG_up",
      "Tanguay_ZF_120hpf_SNOU_up",
      "Tanguay_ZF_120hpf_SOMI_up",
      "Tanguay_ZF_120hpf_SWIM_up",
      "Tanguay_ZF_120hpf_TR_up",
      "Tanguay_ZF_120hpf_TRUN_up",
      "Tanguay_ZF_120hpf_YSE_up"
    )] <- "zebrafish"

  end_point_info$intended_target_family[is.na(end_point_info$intended_target_family)] <- "Undefined"

  cleanUpNames <- end_point_info$intended_target_family
  cleanUpNames <- tools::toTitleCase(cleanUpNames)
  cleanUpNames[grep("Dna", cleanUpNames)] <- "DNA Binding"
  cleanUpNames[grep("Cyp", cleanUpNames)] <- "CYP"
  cleanUpNames[grep("Gpcr", cleanUpNames)] <- "GPCR"
  end_point_info$intended_target_family <- cleanUpNames

  end_point_info$intended_target_family_sub[end_point_info$assay_component_endpoint_name %in%
    c(
      "CLD_CYP1A1_24hr", "CLD_CYP1A1_48hr",
      "CLD_CYP1A1_6hr", "CLD_CYP1A2_24hr",
      "CLD_CYP1A2_48hr", "CLD_CYP1A2_6hr"
    )] <- "basic helix-loop-helix protein"

  end_point_info$intended_target_family_sub[end_point_info$assay_component_endpoint_name %in%
    c(
      "CLD_CYP2B6_24hr", "CLD_CYP2B6_48hr",
      "CLD_CYP2B6_6hr", "CLD_CYP3A4_24hr",
      "CLD_CYP3A4_48hr", "CLD_CYP3A4_6hr",
      "CLD_SULT2A_48hr", "CLD_UGT1A1_48hr"
    )] <- "non-steroidal"

  end_point_info$intended_target_family_sub[end_point_info$assay_component_endpoint_name %in%
    c("NVS_NR_bER", "NVS_NR_bPR", "NVS_NR_cAR")] <- "steroidal"


  cleanUpNames <- end_point_info$intended_target_family_sub
  cleanUpNames <- tools::toTitleCase(cleanUpNames)
  cleanUpNames[grep("Dna", cleanUpNames)] <- "DNA Conformation"
  end_point_info$intended_target_family_sub <- cleanUpNames

  return(end_point_info)
}
