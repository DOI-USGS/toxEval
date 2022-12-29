context("Data")


test_that("Check included data", {
  
  ToxCast_ACC <- ToxCast_ACC
  expect_true(all(names(ToxCast_ACC) %in% c("CAS",
                                            "flags",
                                            "endPoint",
                                            "ACC")))
  
  expect_true(is.numeric(ToxCast_ACC$ACC))
  expect_true(is.character(ToxCast_ACC$CAS))
  expect_true(is.character(ToxCast_ACC$flags))
  expect_true(is.character(ToxCast_ACC$endPoint))
  
  CAS <- unique(ToxCast_ACC$CAS)
  ACC <- get_ACC(CAS)
  
  all_flags <- c(
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
  )
  
  all_gone <- remove_flags(ACC, all_flags)
  expect_true(all(is.na(all_gone$flags)))
  
  end_point_info <- end_point_info
  
  expect_true(all(c("assay_source_name",
                    "assay_component_endpoint_name",
                    "intended_target_family") %in% names(end_point_info)))
  
  # Columns needed for ToxMixtures:
  tm_cols <- c("assay_source_name",
               "assay_name",
               "assay_component_name",
               "assay_component_endpoint_name",
               "intended_target_gene_id",
               "intended_target_gene_name",
               "intended_target_gene_symbol",
               "signal_direction",
               "analysis_direction",
               "biological_process_target")
  
  expect_true(all(tm_cols %in% names(end_point_info)))
  
  default_eps <- c("ACEA", "APR", "ATG",
                   "NVS", "OT", "TOX21", "CEETOX",
                   "LTEA", "CLD", "TANGUAY", "CCTE_PADILLA",
                   "CCTE", "STM", "ARUNA", "CCTE_SHAFER",
                   "CPHEA_STOKER", "CCTE_GLTED", "UPITT", "UKN",
                   "ERF", "TAMU", "IUF", "CCTE_MUNDY", "UTOR", "VALA")
  
  expect_true(all(unique(end_point_info$assay_source_name) %in%
                    c(default_eps, "BSK"))) 
  
  tox_chemicals <- tox_chemicals
  expect_true(all(c("Substance_CASRN",
                    "Structure_MolWt",
                    "Substance_Name",
                    "Total_tested",
                    "Active",
                    "min_concentration",
                    "max_concentration") %in% names(tox_chemicals)))
  
})