#' clean_endPoint_info
#' 
#' Clean up the endPointInfo table from toxCast. Filtering and cleaning based on ES&T (cite Dan/Brett's paper)
#' 
#' @param endPointInfo data frame Endpoint information from ToxCast
#' @export
#' @importFrom stringi stri_trans_totitle
#' @examples 
#' endPointInfo <- endPointInfo
#' cleaned_ep <- clean_endPoint_info(endPointInfo)
clean_endPoint_info <- function(endPointInfo){
  
  endPointInfo <- endPointInfo[!(endPointInfo$assay_source_name == "ATG" & endPointInfo$signal_direction == "loss"),]
  endPointInfo <- endPointInfo[!(endPointInfo$assay_source_name == "NVS" & endPointInfo$signal_direction == "gain"),]
  endPointInfo <- endPointInfo[endPointInfo$assay_component_endpoint_name != "TOX21_p53_BLA_p3_ratio",]
  endPointInfo <- endPointInfo[endPointInfo$assay_component_endpoint_name != "TOX21_p53_BLA_p2_viability",]
  
  endPointInfo$intended_target_family[endPointInfo$assay_component_endpoint_name %in% 
                                        c("CLD_CYP1A1_24hr","CLD_CYP1A1_48hr","CLD_CYP1A1_6hr",
                                          "CLD_CYP1A2_24hr","CLD_CYP1A2_48hr","CLD_CYP1A2_6hr")] <- "dna binding"
  
  endPointInfo$intended_target_family[endPointInfo$assay_component_endpoint_name %in% 
                                        c("CLD_CYP2B6_24hr","CLD_CYP2B6_48hr","CLD_CYP2B6_6hr",
                                          "CLD_CYP3A4_24hr","CLD_CYP3A4_48hr","CLD_CYP3A4_6hr",
                                          "CLD_SULT2A_48hr","CLD_UGT1A1_48hr","NVS_NR_bER",
                                          "NVS_NR_bPR","NVS_NR_cAR")] <- "nuclear receptor"
  
  endPointInfo$intended_target_family[endPointInfo$assay_component_endpoint_name %in% 
                                        c("Tanguay_ZF_120hpf_ActivityScore",
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
                                          "Tanguay_ZF_120hpf_YSE_up")] <- "zebrafish"
  
  endPointInfo$intended_target_family[is.na(endPointInfo$intended_target_family)] <- "Undefined"
  
  assays <- c("ATG","NVS","OT","TOX21","CEETOX", "APR", # only "BSK" is out
              "CLD","TANGUAY","NHEERL_PADILLA",
              "NCCT_SIMMONS","ACEA")
  
  cleanUpNames <- endPointInfo$intended_target_family
  cleanUpNames <- stringi::stri_trans_totitle(cleanUpNames)
  cleanUpNames[grep("Dna",cleanUpNames)] <- "DNA Binding"
  cleanUpNames[grep("Cyp",cleanUpNames)] <- "CYP"
  cleanUpNames[grep("Gpcr",cleanUpNames)] <- "GPCR"
  endPointInfo$intended_target_family <- cleanUpNames

  endPointInfo$intended_target_family_sub[endPointInfo$assay_component_endpoint_name %in% 
                                            c("CLD_CYP1A1_24hr","CLD_CYP1A1_48hr",
                                              "CLD_CYP1A1_6hr","CLD_CYP1A2_24hr",
                                              "CLD_CYP1A2_48hr","CLD_CYP1A2_6hr")] <- "basic helix-loop-helix protein"
  
  endPointInfo$intended_target_family_sub[endPointInfo$assay_component_endpoint_name %in% 
                                            c("CLD_CYP2B6_24hr","CLD_CYP2B6_48hr",
                                              "CLD_CYP2B6_6hr","CLD_CYP3A4_24hr",
                                              "CLD_CYP3A4_48hr","CLD_CYP3A4_6hr",
                                              "CLD_SULT2A_48hr","CLD_UGT1A1_48hr")] <- "non-steroidal"
  
  endPointInfo$intended_target_family_sub[endPointInfo$assay_component_endpoint_name %in% 
                                            c("NVS_NR_bER","NVS_NR_bPR","NVS_NR_cAR")] <- "steroidal"
  
  
  cleanUpNames <- endPointInfo$intended_target_family_sub
  cleanUpNames <- stri_trans_totitle(cleanUpNames)
  cleanUpNames[grep("Dna",cleanUpNames)] <- "DNA Conformation"
  endPointInfo$intended_target_family_sub <- cleanUpNames
  
  return(endPointInfo)
}