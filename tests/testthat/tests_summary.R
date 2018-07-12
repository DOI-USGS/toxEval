context("Chemical summary tests")

path_to_tox <-  system.file("extdata", package="toxEval")
file_name <- "OWC_data_fromSup.xlsx"
full_path <- file.path(path_to_tox, file_name)

chem_data <- readxl::read_excel(full_path, sheet = "Data")
chem_info <- readxl::read_excel(full_path, sheet = "Chemicals")
chem_site <- readxl::read_excel(full_path, sheet = "Sites")
exclusion <- readxl::read_excel(full_path, sheet = "Exclude")

ACClong <- get_ACC(chem_info$CAS)
ACClong <- remove_flags(ACClong)

cleaned_ep <- clean_endPoint_info(endPointInfo)
filtered_ep <- filter_groups(cleaned_ep)

tox_list <- create_toxEval(full_path)

chemicalSummary1 <- get_chemical_summary(tox_list = NULL,
                                         ACClong,
                                         filtered_ep,
                                         chem_data,
                                         chem_site,
                                         chem_info,
                                         exclusion)

chemicalSummary <- get_chemical_summary(tox_list,ACClong,filtered_ep)

test_that("Calculating tox_list", {
  
  expect_type(tox_list, "list")
  expect_length(tox_list, 5)
  expect_equivalent(chemicalSummary, chemicalSummary1)
  
  expect_warning(tox_list_2 <- create_toxEval(file.path(path_to_tox, "test_data.xlsx")))
  expect_true(all(c("SiteID","Sample Date","CAS","Value") %in% names(tox_list_2$chem_data)))
  expect_true(all(c("Class","CAS") %in% names(tox_list_2$chem_info)))
  expect_true(all(c("SiteID","dec_lat","dec_lon","Short Name") %in% names(tox_list_2$chem_site)))
  
  
})

test_that("Calculating summaries", {
  testthat::skip_on_cran()
 
  expect_true(all(names(chemicalSummary) %in% c("CAS","chnm","endPoint","site","date",
                                                "EAR","shortName","Class","Bio_category")))
  
  expect_gte(min(chemicalSummary$EAR), 0)
  expect_true(all(!is.na(chemicalSummary$EAR)))
  
  expect_true(all(chemicalSummary$endPoint %in% ACClong$endPoint))
  
})

test_that("Plotting summaries", {
  testthat::skip_on_cran()
  
  bioPlot <- plot_tox_boxplots(chemicalSummary, 
                               category = "Biological")

  expect_true(all(names(bioPlot$data) %in% c("site","category","meanEAR")))
  expect_equal(bioPlot$layers[[2]]$geom_params$outlier.shape, 19)
  expect_equal(bioPlot$layers[[2]]$aes_params$fill, "steelblue")
  
  classPlot <- plot_tox_boxplots(chemicalSummary, 
                               category = "Chemical Class")

  expect_true(all(names(classPlot$data) %in% c("site","category","meanEAR")))
  expect_equal(classPlot$layers[[2]]$geom_params$outlier.shape, 19)
  expect_equal(classPlot$layers[[2]]$aes_params$fill, "steelblue")
  
  chemPlot <- suppressWarnings(plot_tox_boxplots(chemicalSummary, 
                                 category = "Chemical"))
  
  expect_true(all(names(chemPlot$data) %in% c("site","chnm","Class","meanEAR")))
  expect_equal(chemPlot$layers[[1]]$geom_params$outlier.shape, 19)
  expect_true(is.null(chemPlot$layers[[1]]$aes_params$fill))
  
})

test_that("Plotting heatmap summaries", {
  testthat::skip_on_cran()
  
  bioHeatPlot <- plot_tox_heatmap(chemicalSummary, chem_site,
                               category = "Biological")
  expect_true(all(names(bioHeatPlot$data) %in% c("site","category","meanEAR",
                                                 "site_grouping","Short Name")))

  classHeatPlot <- plot_tox_heatmap(chemicalSummary, chem_site,
                                  category = "Chemical Class")
  expect_true(all(names(classHeatPlot$data) %in% c("site","category","meanEAR",
                                                 "site_grouping","Short Name")))
  
  chemHeatPlot <- plot_tox_heatmap(chemicalSummary, chem_site,
                                    category = "Chemical")
  expect_true(all(names(chemHeatPlot$data) %in% c("site","chnm","Class","meanEAR",
                                                   "site_grouping","Short Name")))
  
})

test_that("Plotting stacked summaries", {
  testthat::skip_on_cran()
  
  bioStackPlot <- plot_tox_stacks(chemicalSummary, chem_site,
                                  category = "Biological")
  expect_true(all(names(bioStackPlot$data) %in% c("site","category","meanEAR",
                                                 "site_grouping","Short Name")))
  
  classStackPlot <- plot_tox_stacks(chemicalSummary, chem_site,
                                    category = "Chemical Class")
  expect_true(all(names(classStackPlot$data) %in% c("site","category","meanEAR",
                                                   "site_grouping","Short Name")))
  
  chemStackPlot <- plot_tox_stacks(chemicalSummary, chem_site,
                                   category = "Chemical",include_legend = FALSE)
  expect_true(all(names(chemStackPlot$data) %in% c("site","category","meanEAR",
                                                  "site_grouping","Short Name","Class")))
  
})

test_that("Plotting endpoints", {
  testthat::skip_on_cran()
  
  bioStackPlot <- suppressWarnings(plot_tox_endpoints(chemicalSummary, 
                                  category = "Biological",
                                  filterBy = "Cell Cycle"))
  expect_true(all(names(bioStackPlot$data) %in% c("site","category","endPoint",
                                                  "meanEAR")))
  
  classStackPlot <- suppressWarnings(plot_tox_endpoints(chemicalSummary, 
                                    category = "Chemical Class", filterBy = "PAHs"))
  expect_true(all(names(classStackPlot$data) %in% c("site","category","endPoint",
                                                    "meanEAR")))
  
  chemStackPlot <- suppressWarnings(plot_tox_endpoints(chemicalSummary, category = "Chemical", filterBy = "Atrazine"))
  expect_true(all(names(chemStackPlot$data) %in% c("site","category","endPoint",
                                                   "meanEAR")))
  
})


test_that("Table functions", {
  testthat::skip_on_cran()
  
  statStuff <- rank_sites(chemicalSummary, "Biological", hit_threshold = 0.1, mean_logic = FALSE)

  expect_true(all(c("site","DNA Binding maxEAR",     
                    "DNA Binding freq","Nuclear Receptor maxEAR",
                    "Nuclear Receptor freq","Esterase maxEAR",
                    "Cell Cycle freq", "Cell Cycle maxEAR") %in% names(statStuff)))
  expect_equal(round(statStuff[["DNA Binding maxEAR"]][which(statStuff[["site"]] == "Raisin")],4),4.2992)
  expect_equal(round(statStuff[["DNA Binding freq"]][which(statStuff[["site"]] == "Raisin")],4),0.8409)
  
  groupStuff <- hits_summary(chemicalSummary, "Biological", hit_threshold = 1)
  expect_true(all(unique(groupStuff$site) %in% chem_site$`Short Name`))
  expect_true(all(c("site","category","Samples with hits","Number of Samples") %in% names(groupStuff)))
  expect_equal(groupStuff[["Samples with hits"]][which(groupStuff[["site"]] == "Raisin" &
                                                           groupStuff[["category"]] == "DNA Binding")],28)
  expect_equal(groupStuff[["Number of Samples"]][which(groupStuff[["site"]] == "Raisin" &
                                                       groupStuff[["category"]] == "DNA Binding")],44)
  
})

test_that("Chem plotting functions", {
  testthat::skip_on_cran()
  
  graphData <- graph_chem_data(chemicalSummary)
  expect_true(all(names(graphData) %in% c("site","chnm","Class","meanEAR")))
  expect_equal(levels(graphData$Class)[1], "Detergent Metabolites")
  expect_equal(levels(graphData$Class)[length(levels(graphData$Class))], "Fuels")
  expect_equal(signif(graphData[["meanEAR"]][graphData[["site"]] == "USGS-04024000" &
                                       graphData[["chnm"]] == "1-Methylnaphthalene"],4),2.279e-06)
  expect_equal(signif(graphData[["meanEAR"]][graphData[["site"]] == "USGS-04024000" &
                                              graphData[["chnm"]] == "4-Nonylphenol, branched"],4),0.9975)
  
  
})

test_that("Map stuff functions", {
  testthat::skip_on_cran()
  
  mapDataList <- map_tox_data(chemicalSummary, 
                            chem_site = chem_site, 
                            category = "Biological")
  expect_type(mapDataList, "list")
  expect_equal(length(mapDataList), 2)
  map_df <- mapDataList[["mapData"]]
  expect_equal(signif(map_df[["meanMax"]][map_df[["Short Name"]] == "StLouis"],4),1.271)
  expect_equal(map_df[["count"]][map_df[["Short Name"]] == "StLouis"],26)
  expect_equal(map_df[["sizes"]][map_df[["Short Name"]] == "StLouis"],6.6)
  
  # Map data:
  map_data <- make_tox_map(chemicalSummary, chem_site, "Biological")
  expect_type(map_data, "list")
  expect_equal(length(map_data), 8) 
  expect_true(all(class(map_data) %in% c("leaflet","htmlwidget")))
})

test_that("Table endpoint hits", {
  testthat::skip_on_cran()

  bt <- endpoint_hits_DT(chemicalSummary, category = "Biological")
  expect_type(bt, "list")
  expect_true(all(names(bt$x$data) %in% c("endPoint","Nuclear Receptor","DNA Binding",
                                  "Phosphatase","Steroid Hormone","Esterase")))
  expect_true(all(class(bt) %in% c("datatables","htmlwidget")))
  
  bt_df <- endpoint_hits(chemicalSummary, category = "Biological")
  expect_true(all(names(bt_df) %in% c("endPoint","Nuclear Receptor","DNA Binding",     
                                      "Esterase","Steroid Hormone")))
  expect_equal(bt_df[["Nuclear Receptor"]][bt_df[["endPoint"]] == "NVS_NR_hPPARg"],11)
  expect_true(is.na(bt_df[["Esterase"]][bt_df[["endPoint"]] == "NVS_NR_hPPARg"]))
  
  expect_error(endpoint_hits_DT(chemicalSummary, category = "Class"))
  
  ct <- endpoint_hits_DT(chemicalSummary, category = "Chemical Class")
  expect_type(ct, "list")
  
  expect_true(all(names(ct$x$data) %in% c("endPoint","Antioxidants","PAHs",
                                          "Detergent Metabolites","Herbicides",
                                          "Plasticizers")))
  
  cht <- endpoint_hits_DT(chemicalSummary, category = "Chemical")
  expect_type(cht, "list")
  
  expect_true(all(names(cht$x$data) %in% c("endPoint","Bisphenol A","Fluoranthene","Tris(2-chloroethyl) phosphate",
                                           "4-Nonylphenol, branched","Triphenyl phosphate",
                                           "Metolachlor","Atrazine")))
})

test_that("hits_by_groupings_DT", {
  testthat::skip_on_cran()
  
  bt <- hits_by_groupings_DT(chemicalSummary, category = "Biological")
  expect_type(bt, "list")
  expect_true(all(class(bt) %in% c("datatables","htmlwidget")))
  expect_true("nSites" %in% names(bt$x$data))
  
  bt_df <- hits_by_groupings(chemicalSummary, category = "Chemical Class")
  expect_true(all(names(bt_df) %in% c("Nuclear Receptor","DNA Binding","Esterase",        
                                      "Steroid Hormone","Zebrafish")))
  expect_true(all(c("Detergent Metabolites","Antioxidants","Herbicides") %in% rownames(bt_df)))
  expect_equal(bt_df[["Nuclear Receptor"]], c(10,11,7,11,rep(0,9)))
  
  expect_error(hits_by_groupings_DT(chemicalSummary, category = "Class"))
  
  ct <- hits_by_groupings_DT(chemicalSummary, category = "Chemical Class")
  expect_type(ct, "list")
  expect_true(all(class(ct) %in% c("datatables","htmlwidget")))
  expect_true(all(names(ct$x$data) %in% c(" ","Nuclear Receptor","DNA Binding",
                                          "Phosphatase","Esterase","Steroid Hormone",
                                          "Zebrafish")))
  
  cht <- hits_by_groupings_DT(chemicalSummary, category = "Chemical")
  expect_type(cht, "list")
  
  expect_true(all(names(cht$x$data) %in% c(" ","Nuclear Receptor","DNA Binding",
                                           "Phosphatase","Esterase","Steroid Hormone",
                                           "Zebrafish")))
})


test_that("hits_summary_DT", {
  testthat::skip_on_cran()
  
  bt <- hits_summary_DT(chemicalSummary, category = "Biological")
  expect_type(bt, "list")
  expect_true(all(class(bt) %in% c("datatables","htmlwidget")))
  expect_true(all(c("site","category","Samples with hits","Number of Samples") %in% names(bt$x$data)))
  
  expect_error(hits_summary_DT(chemicalSummary, category = "Class"))
  
  ct <- hits_summary_DT(chemicalSummary, category = "Chemical Class")
  expect_type(ct, "list")
  
  expect_true(all(names(ct$x$data) %in% c("site","category","Samples with hits","Number of Samples")))
  
  cht <- hits_summary_DT(chemicalSummary, category = "Chemical")
  expect_type(cht, "list")
  
  expect_true(all(names(cht$x$data) %in% c("site","category","Samples with hits","Number of Samples")))
})

test_that("rank_sites_DT", {
  testthat::skip_on_cran()
  
  bt <- rank_sites_DT(chemicalSummary, category = "Biological")
  expect_type(bt, "list")
  expect_true("site" %in% names(bt$x$data))
  
  expect_error(rank_sites_DT(chemicalSummary, category = "Class"))
  
  ct <- rank_sites_DT(chemicalSummary, category = "Chemical Class")
  expect_type(ct, "list")
  
  expect_true(all(c("site","Antioxidants maxEAR","Antioxidants freq",
                    "Herbicides maxEAR","Herbicides freq",
                    "Detergent Metabolites maxEAR","Detergent Metabolites freq") %in% names(ct$x$data)))
  
  cht <- rank_sites_DT(chemicalSummary, category = "Chemical")
  expect_type(cht, "list")
  
  expect_true(all(c("site","Bisphenol A maxEAR","Bisphenol A freq",
                    "4-Nonylphenol, branched maxEAR") %in% names(cht$x$data)))
})