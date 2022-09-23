
.onAttach <- function(libname, pkgname) {

  packageStartupMessage(
    paste(strwrap(paste('USGS Research Package:
https://owi.usgs.gov/R/packages.html#research
ToxCast database: version', dbVersion()), width = 40),
      collapse='\n'))
}

dbVersion <- function(){
  "3.5"
}

#' Analyze ToxCast data in relation to measured concentrations.
#' 
#' \code{toxEval} includes a set of functions to analyze, visualize, and 
#' organize measured concentration data as it relates to ToxCast data 
#' (default) or other user-selected chemical-biological interaction 
#' benchmark data such as water quality criteria. The intent of 
#' these analyses is to develop a better understanding of the potential 
#' biological relevance of environmental chemistry data. Results can 
#' be used to prioritize which chemicals at which sites may be of 
#' greatest concern. These methods are meant to be used as a screening 
#' technique to predict potential for biological influence from chemicals 
#' that ultimately need to be validated with direct biological assays. 

#'
#' \tabular{ll}{
#' Package: \tab toxEval\cr
#' Type: \tab Package\cr
#' License: \tab Unlimited for this package, dependencies have more restrictive licensing.\cr
#' Copyright: \tab This software is in the public domain because it contains materials
#' that originally came from the United States Geological Survey, an agency of
#' the United States Department of Interior. For more information, see the
#' official USGS copyright policy at
#' https://www.usgs.gov/visual-id/credit_usgs.html#copyright\cr
#' LazyLoad: \tab yes\cr
#' }
#'
#'
#' @name toxEval-package
#' @docType package
#' @author Laura De Cicco \email{ldecicco@@usgs.gov}. Steven Corsi  
#' @keywords ToxCast
NULL

#' ACC values included with toxEval. 
#' 
#' Downloaded on October 2020 from ToxCast. The data were
#' combined from files in the "INVITRODB_V3_3_LEVEL5" folder. 
#' At the time of toxEval package release, this information was found:
#' \url{https://www.epa.gov/chemical-research/exploring-toxcast-data-downloadable-data}
#' in the "ToxCast & Tox21 Data Spreadsheet" data set. 
#' ACC values are the in the "ACC" column (winning model) and units are 
#' log micro-Molarity (log \eqn{\mu}M).
#' 
#' @references Toxicology, EPA's National Center for Computational (2020): ToxCast and Tox21 Data Spreadsheet. figshare. Dataset.
#'  https://doi.org/10.23645/epacomptox.6062503.v3
#'  
#' @source \url{https://www.epa.gov/chemical-research/exploring-toxcast-data-downloadable-data}
#'
#'@aliases ToxCast_ACC
#'@return data frame with columns CAS, chnm (chemical name), flags, endPoint, and ACC (value).
#'@name ToxCast_ACC
#'@docType data
#'@export ToxCast_ACC
#'@keywords datasets
#'@examples
#'head(ToxCast_ACC)
NULL


# This is staying commented out because it adds 
# extraneous notes:
# path_to_files <- "../toxCast_Data/INVITRODB_V3_3_LEVEL5"

# files <- list.files(path = path_to_files)
# 
# x <- data.table::fread(file.path(path_to_files, files[1]), data.table = FALSE)
# 
# filtered <- dplyr::select(x, chnm, casn, aenm, logc_min, logc_max, modl_acc,
#                           modl, actp, modl_ga, flags, hitc,gsid_rep)
# filtered <- dplyr::filter(x, hitc == 1)
# 
# for(i in files[-1]){
#   subX <- data.table::fread(file.path(path_to_files,i), data.table = FALSE)
# 
#   subFiltered <- dplyr::select(subX, chnm, casn, aenm, logc_min, logc_max, modl_acc,
#                                modl, actp, modl_ga, flags, hitc,gsid_rep)
#   subFiltered <- dplyr::filter(subFiltered, hitc == 1)
# 
#   filtered <- dplyr::bind_rows(filtered, subFiltered)
# }
# 
# ACCgain <- dplyr::filter(filtered, hitc == 1)
# ACCgain <- dplyr::filter(ACCgain, gsid_rep == 1)
# ACCgain <- dplyr::select(ACCgain, casn, chnm, aenm, modl_acc, flags)
# ACCgain <- tidyr::spread(ACCgain, key = aenm, value = modl_acc)
# 
# ACC <- ACCgain
# ACC <- tidyr::gather(ACC, endPoint, ACC, -casn, -chnm, -flags)

# ACC <- filter(ACC, !is.na(ACC))
# ACC <- rename(ACC, CAS = casn)
# ACC <- select(ACC, -chnm)
# 
# saveRDS(ACC, "ACC_v3.rds")
# ToxCast_ACCv3 <- readRDS("ACC_v3.rds")
# 
# 
# end_point_info_v3_assay <- readxl::read_xlsx("../toxCast_Data/INVITRODB_V3_3_SUMMARY/assay_annotation_information_invitrodb_v3_3.xlsx", sheet = "assay")
# 
# end_point_info_v3_assay.component <- readxl::read_xlsx("../toxCast_Data/INVITRODB_V3_3_SUMMARY/assay_annotation_information_invitrodb_v3_3.xlsx",
#                                                        sheet = "assay.component")
# end_point_info_v3_assay.component.endpoint <- readxl::read_xlsx("../toxCast_Data/INVITRODB_V3_3_SUMMARY/assay_annotation_information_invitrodb_v3_3.xlsx",
#                                                                 sheet = "assay.component.endpoint")
# 
# library(dplyr)
# 
# end_point_info_v3 <- end_point_info_v3_assay %>%
#   left_join(end_point_info_v3_assay.component, by = "aid") %>%
#   left_join(end_point_info_v3_assay.component.endpoint, by = "acid")
# 
# gene_stuff <- readxl::read_xlsx("../toxCast_Data/INVITRODB_V3_3_SUMMARY/gene_target_information_invitrodb_v3_3.xlsx")
# 
#Generate table for EPA review. Look at endpoint-gene linkages and assess validity
# gene_stuff_out <- gene_stuff %>%
#   select(gene_symbol, gene_name,
#          aeid, aenm) %>%
#   distinct() %>%
#   arrange(gene_symbol, aenm) %>%
#   left_join(end_point_info_v3, by = c("aeid", "aenm" = "assay_component_endpoint_name")) %>%
#   select(gene_symbol, gene_name,
#          aeid, aenm, signal_direction, analysis_direction, assay_component_endpoint_desc)
# 
# write.csv(gene_stuff_out, "genes_endpoints_for_review.csv", row.names = FALSE)
# After reviewing the gene-endpoint combos, some additional changes will need to happen to the gene table
# removing some gene-endpoint linkages. Add additional columns to note agonism/antagomism. Etc.

# Merge gene with rest of endpoint table,
# Collapse gene table to "one row = one endpoint"
# endpoints with multiple genes will have |'s in gene columns.
# 
# gene_stuff2 <- gene_stuff %>%
#   select(-organism_id) %>%
#   group_by(aeid, aenm) %>%
#   summarize(across(everything(), function(x) paste(unique(x[which(x != "" & x != "NA" & !is.na(x) )]),
#                                                    collapse = "|")), .groups = "drop") %>%
#   rename(intended_target_gene_id = gene_id,
#          intended_target_gene_name = gene_name,
#          intended_target_gene_symbol = gene_symbol)
# 
# end_point_info_v3 <- end_point_info_v3 %>%
#   left_join(gene_stuff2, by = c("aeid", "assay_component_endpoint_name" = "aenm"))
# 
# assay_table <- unique(end_point_info[c("assay_source_name", "assay_source_long_name")])
# 
# end_point_info_v3$assay_source_name <- gsub("\\_.*" , "", end_point_info_v3$assay_name)
# 
# end_point_info_v3$assay_source_name[grepl("NHEERL", end_point_info_v3$assay_source_name)] <-
#   paste0("NHEERL_", (gsub("\\_.*" , "",
#        gsub("NHEERL_", "", end_point_info_v3$assay_name[grepl("NHEERL", end_point_info_v3$assay_source_name)]))))
# 
# end_point_info <- end_point_info_v3
# rm(end_point_info_v3, end_point_info_v3_assay, end_point_info_v3_assay.component,
#    end_point_info_v3_assay.component.endpoint, gene_stuff, gene_stuff2)
#   


#' Endpoint information from ToxCast
#' 
#' Downloaded on October 2020 from ToxCast. The file name of the
#' raw data was "assay_annotation_information_invitrodb_v3_3.xlsx" from the zip file 
#' "INVITRODB_V3_3_SUMMARY" folder. At the time
#' of the toxEval package release, these data were found at:
#' \url{https://www.epa.gov/chemical-research/exploring-toxcast-data-downloadable-data}
#' in the section marked "Download Assay Information", in the 
#' ToxCast & Tox21 high-throughput assay information data set.
#'
#'
#'@name end_point_info
#'@aliases end_point_info
#'@docType data
#'@keywords datasets
#'@references U.S. EPA. 2014. ToxCast Assay Annotation Data User Guide. 
#'\url{https://www.epa.gov/chemical-research/toxcast-assay-annotation-data-user-guide}.
#'  
#'@source \url{https://doi.org/10.23645/epacomptox.6062479.v3}
#'@export end_point_info
#'@return data frame with 86 columns. The columns and definitions
#'are discussed in the "ToxCast Assay Annotation Version 1.0 Data User Guide (PDF)" (see source)
#'@examples
#'end_point_info <- end_point_info
#'head(end_point_info[,1:5])
NULL

#' ToxCast Chemical Information 
#' 
#' Downloaded on October 2015 from ToxCast. The file name of the
#' raw data was "TOX21IDs_v4b_23Oct2014_QCdetails.xlsx", 
#' from the US EPA DSSTox DATA RELEASE OCTOBER 2015. At the time
#' of toxEval package release, this information was found:
#' \url{https://www.epa.gov/chemical-research/exploring-toxcast-data-downloadable-data}
#' in the section marked "Download ToxCast Chemical Information". This 
#' was in the "ToxCast & Tox21 Chemicals Distributed Structure-Searchable Toxicity Database (DSSTox files)"
#' data set.
#'
#'@aliases tox_chemicals
#'@name tox_chemicals
#'@return data frame with columns:    
#'"Substance_Name","Substance_CASRN",
#'"Structure_MolWt" 
#'@docType data
#'@keywords datasets
#'@export tox_chemicals
#'@examples
#'head(tox_chemicals)
NULL
