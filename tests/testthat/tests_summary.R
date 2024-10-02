context("Chemical summary tests")

path_to_tox <- system.file("extdata", package = "toxEval")
file_name <- "OWC_data_fromSup.xlsx"
full_path <- file.path(path_to_tox, file_name)

chem_data <- readxl::read_excel(full_path, sheet = "Data")
chem_info <- readxl::read_excel(full_path, sheet = "Chemicals")
chem_site <- readxl::read_excel(full_path, sheet = "Sites")
exclusion <- readxl::read_excel(full_path, sheet = "Exclude")

ACC <- get_ACC(chem_info$CAS)
ACC <- remove_flags(ACC)

cleaned_ep <- clean_endPoint_info(end_point_info)
filtered_ep <- filter_groups(cleaned_ep)

tox_list <- create_toxEval(full_path)

chemical_summary1 <- get_chemical_summary(
  tox_list = NULL,
  ACC,
  filtered_ep,
  chem_data,
  chem_site,
  chem_info,
  exclusion
)

chemical_summary <- get_chemical_summary(tox_list, ACC, filtered_ep)

test_that("Calculating tox_list", {
  expect_type(tox_list, "list")
  expect_length(tox_list, 5)
  expect_equivalent(chemical_summary, chemical_summary1)

  expect_true(all(c("SiteID", "Sample Date", "CAS", "Value") %in% names(tox_list$chem_data)))
  expect_true(all(c("Class", "CAS") %in% names(tox_list$chem_info)))
  expect_true(all(c("SiteID", "dec_lat", "dec_lon", "Short Name") %in% names(tox_list$chem_site)))
})

test_that("Calculating summaries", {
  testthat::skip_on_cran()

  expect_true(all(names(chemical_summary) %in% c(
    "CAS", "chnm", "endPoint", "site", "date",
    "EAR", "shortName", "Class", "Bio_category"
  )))

  expect_gte(min(chemical_summary$EAR), 0)
  expect_true(all(!is.na(chemical_summary$EAR)))

  expect_true(all(chemical_summary$endPoint %in% ACC$endPoint))
})

test_that("Plotting summaries", {
  testthat::skip_on_cran()

  bioPlot <- plot_tox_boxplots(chemical_summary,
    category = "Biological"
  )

  expect_true(all(names(bioPlot$data) %in% c("site", "category", "meanEAR")))
  expect_equal(bioPlot$layers[[2]]$geom_params$outlier.shape, 19)
  expect_equal(bioPlot$layers[[2]]$aes_params$fill, "steelblue")

  classPlot <- plot_tox_boxplots(chemical_summary,
    category = "Chemical Class"
  )

  expect_true(all(names(classPlot$data) %in% c("site", "category", "meanEAR")))
  expect_equal(classPlot$layers[[2]]$geom_params$outlier.shape, 19)
  expect_equal(classPlot$layers[[2]]$aes_params$fill, "steelblue")

  chemPlot <- suppressWarnings(plot_tox_boxplots(chemical_summary,
    category = "Chemical"
  ))

  expect_true(all(names(chemPlot$data) %in% c("site", "chnm", "Class", "meanEAR")))
  expect_equal(chemPlot$layers[[1]]$geom_params$outlier.shape, 19)
  expect_true(is.null(chemPlot$layers[[1]]$aes_params$fill))
})

test_that("Plotting heatmap summaries", {
  testthat::skip_on_cran()

  bioHeatPlot <- plot_tox_heatmap(chemical_summary, chem_site,
    category = "Biological"
  )
  expect_true(all(names(bioHeatPlot$data) %in% c(
    "site", "category", "meanEAR",
    "site_grouping", "Short Name"
  )))

  classHeatPlot <- plot_tox_heatmap(chemical_summary, chem_site,
    category = "Chemical Class"
  )
  expect_true(all(names(classHeatPlot$data) %in% c(
    "site", "category", "meanEAR",
    "site_grouping", "Short Name"
  )))

  chemHeatPlot <- plot_tox_heatmap(chemical_summary, chem_site,
    category = "Chemical"
  )
  expect_true(all(names(chemHeatPlot$data) %in% c(
    "site", "chnm", "Class", "meanEAR",
    "site_grouping", "Short Name"
  )))
})

test_that("Plotting stacked summaries", {
  testthat::skip_on_cran()

  bioStackPlot <- plot_tox_stacks(chemical_summary, chem_site,
    category = "Biological"
  )
  expect_true(all(names(bioStackPlot$data) %in% c(
    "site", "category", "meanEAR",
    "site_grouping", "Short Name"
  )))

  classStackPlot <- plot_tox_stacks(chemical_summary, chem_site,
    category = "Chemical Class"
  )
  expect_true(all(names(classStackPlot$data) %in% c(
    "site", "category", "meanEAR",
    "site_grouping", "Short Name"
  )))

  chemStackPlot <- plot_tox_stacks(chemical_summary, chem_site,
    category = "Chemical", include_legend = FALSE
  )
  expect_true(all(names(chemStackPlot$data) %in% c(
    "site", "category", "meanEAR",
    "site_grouping", "Short Name", "Class"
  )))
})

test_that("Plotting endpoints", {
  testthat::skip_on_cran()

  bioStackPlot <- suppressWarnings(plot_tox_endpoints(chemical_summary,
    category = "Biological",
    filterBy = "Cell Cycle"
  ))
  expect_true(all(names(bioStackPlot$data) %in% c(
    "site", "category", "endPoint",
    "meanEAR"
  )))

  classStackPlot <- suppressWarnings(plot_tox_endpoints(chemical_summary,
    category = "Chemical Class", filterBy = "PAHs"
  ))
  expect_true(all(names(classStackPlot$data) %in% c(
    "site", "category", "endPoint",
    "meanEAR"
  )))

  chemStackPlot <- suppressWarnings(plot_tox_endpoints(chemical_summary, category = "Chemical", filterBy = "Atrazine"))
  expect_true(all(names(chemStackPlot$data) %in% c(
    "site", "category", "endPoint",
    "meanEAR"
  )))
})


test_that("Table functions", {
  testthat::skip_on_cran()

  statStuff <- rank_sites(chemical_summary, "Biological", hit_threshold = 0.1, mean_logic = FALSE)

  expect_true(all(c(
    "site", "DNA Binding maxEAR",
    "DNA Binding freq", "Nuclear Receptor maxEAR",
    "Nuclear Receptor freq", "Esterase maxEAR",
    "Cell Cycle freq", "Cell Cycle maxEAR"
  ) %in% names(statStuff)))

  expect_equal(round(statStuff[["DNA Binding maxEAR"]][which(statStuff[["site"]] == "Raisin")], 4), 0.0301)

  expect_equal(round(statStuff[["DNA Binding freq"]][which(statStuff[["site"]] == "Raisin")], 4), 0)

  groupStuff <- hits_summary(chemical_summary, "Biological", hit_threshold = 1)
  expect_true(all(unique(groupStuff$site) %in% chem_site$`Short Name`))
  expect_true(all(c("site", "category", "Samples with hits", "Number of Samples") %in% names(groupStuff)))
  expect_equal(groupStuff[["Samples with hits"]][which(groupStuff[["site"]] == "Raisin" &
    groupStuff[["category"]] == "Nuclear Receptor")], 29)

  expect_equal(groupStuff[["Number of Samples"]][which(groupStuff[["site"]] == "Raisin" &
    groupStuff[["category"]] == "DNA Binding")], 44)
})

test_that("Chem plotting functions", {
  testthat::skip_on_cran()

  graphData <- graph_chem_data(chemical_summary)
  expect_true(all(names(graphData) %in% c("site", "chnm", "Class", "meanEAR")))
  expect_equal(levels(graphData$Class)[1], "Herbicides")
  expect_equal(levels(graphData$Class)[length(levels(graphData$Class))], "Dyes/Pigments")
  expect_equal(signif(graphData[["meanEAR"]][graphData[["site"]] == "USGS-04024000" &
    graphData[["chnm"]] == "Naphthalene"], 4), 0.0004799)
  expect_equal(signif(graphData[["meanEAR"]][graphData[["site"]] == "USGS-04024000" &
    graphData[["chnm"]] == "Anthraquinone"], 4), 6.131e-05)
})

test_that("Map stuff functions", {
  testthat::skip_on_cran()

  mapDataList <- map_tox_data(chemical_summary,
    chem_site = chem_site,
    category = "Biological"
  )
  expect_type(mapDataList, "list")
  expect_equal(length(mapDataList), 2)
  map_df <- mapDataList[["mapData"]]
  expect_equal(signif(map_df[["meanMax"]][map_df[["Short Name"]] == "StLouis"], 4), 1.891)

  expect_equal(map_df[["count"]][map_df[["Short Name"]] == "StLouis"], 31)
  expect_equal(map_df[["sizes"]][map_df[["Short Name"]] == "StLouis"], 7.2)

  # Map data:
  map_data <- make_tox_map(chemical_summary, chem_site, "Biological")
  expect_type(map_data, "list")
  expect_equal(length(map_data), 8)
  expect_true(all(class(map_data) %in% c("leaflet", "htmlwidget")))
})

test_that("Table endpoint hits", {
  testthat::skip_on_cran()

  bt <- endpoint_hits_DT(chemical_summary, category = "Biological")
  expect_type(bt, "list")
  expect_true(all(names(bt$x$data) %in% c(
    "endPoint", "Nuclear Receptor", "CYP", "Lyase",
    "Steroid Hormone", "Esterase", "Metabolite"
  )))
  expect_true(all(class(bt) %in% c("datatables", "htmlwidget")))

  bt_df <- endpoint_hits(chemical_summary, category = "Biological")
  expect_true(all(names(bt_df) %in% c(
    "endPoint", "Nuclear Receptor", "CYP", "Lyase",
    "Steroid Hormone", "Esterase", "Metabolite"
  )))
  expect_equal(bt_df[["Nuclear Receptor"]][bt_df[["endPoint"]] == "OT_ER_ERaERb_1440"], 1)
  expect_true(is.na(bt_df[["Esterase"]][bt_df[["endPoint"]] == "OT_ER_ERbERb_0480"]))

  expect_error(endpoint_hits_DT(chemical_summary, category = "Class"))

  ct <- endpoint_hits_DT(chemical_summary, category = "Chemical Class")
  expect_type(ct, "list")

  expect_true(all(names(ct$x$data) %in% c(
    "endPoint", "Antioxidants", 
    "Detergent Metabolites", "Herbicides",
    "Antimicrobial Disinfectants"
  )))

  cht <- endpoint_hits_DT(chemical_summary, category = "Chemical")
  expect_type(cht, "list")

  expect_true(all(names(cht$x$data) %in% c(
    "endPoint", "Bisphenol A", "Metolachlor",
    "4-Nonylphenol, branched", "Pyrene",
    "Triclosan", "Atrazine"
  )))
})

test_that("hits_by_groupings_DT", {
  testthat::skip_on_cran()

  bt <- hits_by_groupings_DT(chemical_summary, category = "Biological")
  expect_type(bt, "list")
  expect_true(all(class(bt) %in% c("datatables", "htmlwidget")))
  expect_true("nSites" %in% names(bt$x$data))

  bt_df <- hits_by_groupings(chemical_summary, category = "Chemical Class")
  expect_true(all(names(bt_df) %in% c(
    "Nuclear Receptor", "DNA Binding", "Cell Cycle", "Esterase",
    "Steroid Hormone", "Channel 2", "CYP", "Transporter", "Lyase",
    "Metabolite"
  )))

  expect_true(all(c("Detergent Metabolites", "Antioxidants", "Herbicides") %in% rownames(bt_df)))
  expect_equal(bt_df[["Nuclear Receptor"]], c(8, 10, 16, rep(0, 9), NA))

  expect_error(hits_by_groupings_DT(chemical_summary, category = "Class"))

  ct <- hits_by_groupings_DT(chemical_summary, category = "Chemical Class")
  expect_type(ct, "list")
  expect_true(all(class(ct) %in% c("datatables", "htmlwidget")))
  expect_true(all(names(ct$x$data) %in% c(
    " ", "Nuclear Receptor", "DNA Binding", "Cell Cycle",
    "Esterase", "CYP", "Steroid Hormone", "Lyase",
    "Transporter", "Channel 2", "Metabolite"
  )))


  cht <- hits_by_groupings_DT(chemical_summary, category = "Chemical")
  expect_type(cht, "list")

  expect_true(all(names(cht$x$data) %in% c(
    " ", "Nuclear Receptor", "DNA Binding", "Cell Cycle",
    "CYP", "Esterase", "Steroid Hormone",
    "Channel 2", "Transporter", "Lyase", "Metabolite"
  )))
})


test_that("hits_summary_DT", {
  testthat::skip_on_cran()

  bt <- hits_summary_DT(chemical_summary, category = "Biological")
  expect_type(bt, "list")
  expect_true(all(class(bt) %in% c("datatables", "htmlwidget")))
  expect_true(all(c("site", "category", "Samples with hits", "Number of Samples") %in% names(bt$x$data)))

  expect_error(hits_summary_DT(chemical_summary, category = "Class"))

  ct <- hits_summary_DT(chemical_summary, category = "Chemical Class")
  expect_type(ct, "list")

  expect_true(all(names(ct$x$data) %in% c("site", "category", "Samples with hits", "Number of Samples")))

  cht <- hits_summary_DT(chemical_summary, category = "Chemical")
  expect_type(cht, "list")

  expect_true(all(names(cht$x$data) %in% c("site", "category", "Samples with hits", "Number of Samples")))
})

test_that("rank_sites_DT", {
  testthat::skip_on_cran()

  bt <- rank_sites_DT(chemical_summary, category = "Biological")
  expect_type(bt, "list")
  expect_true("site" %in% names(bt$x$data))

  expect_error(rank_sites_DT(chemical_summary, category = "Class"))

  ct <- rank_sites_DT(chemical_summary, category = "Chemical Class")
  expect_type(ct, "list")

  expect_true(all(c(
    "site", "Antioxidants maxEAR", "Antioxidants freq",
    "Herbicides maxEAR", "Herbicides freq",
    "Detergent Metabolites maxEAR", "Detergent Metabolites freq"
  ) %in% names(ct$x$data)))

  cht <- rank_sites_DT(chemical_summary, category = "Chemical")
  expect_type(cht, "list")

  expect_true(all(c(
    "site", "Bisphenol A maxEAR", "Bisphenol A freq",
    "4-Nonylphenol, branched maxEAR"
  ) %in% names(cht$x$data)))
})

test_that("Calculating completness", {
  testthat::skip_on_cran()

  graphData <- graph_chem_data(chemical_summary)
  complete_df <- toxEval:::get_complete_set(chemical_summary, graphData, tox_list$chem_site)
  expect_equal(length(unique(complete_df$chnm)), length(levels(chemical_summary$chnm)))
  expect_equal(nrow(complete_df), length(levels(chemical_summary$chnm)) * length(unique(tox_list$chem_site$`Short Name`)))

  graphData2 <- tox_boxplot_data(chemical_summary, "Biological")
  complete_df_cat <- toxEval:::get_complete_set_category(chemical_summary, graphData2, tox_list$chem_site, category = "Biological")

  expect_equal(nrow(complete_df_cat), 2451)
})

test_that("Calculating concentrations", {
  testthat::skip_on_cran()

  summary_conc <- get_concentration_summary(tox_list)

  expect_equal(nrow(tox_list$chem_data), nrow(summary_conc))
  expect_equivalent(summary_conc$EAR, chem_data$Value)

  gd_conc <- graph_chem_data(summary_conc)
  gd_tox <- graph_chem_data(chemical_summary)

  combo <- side_by_side_data(gd_conc, gd_tox,
    left_title = "Concentration",
    right_title = "ToxCast"
  )

  expect_true("guide_side" %in% names(combo))
  expect_equal(
    c("Concentration", "ToxCast"),
    levels(combo$guide_side)
  )
})

test_that("Testing levels", {
  chem_levels <- levels(chemical_summary$chnm)
  class_levels <- levels(chemical_summary$Class)

  gd <- graph_chem_data(chemical_summary)

  expect_equal(levels(gd$chnm), chem_levels)
  expect_equal(levels(gd$Class), class_levels)


  expect_equal(chem_levels[1:5], c(
    "Anthraquinone", "Tetrachloroethylene", "Isophorone",         
    "1,4-Dichlorobenzene", "Bromoform"
  ))
  expect_equal(class_levels[1:5], c(
    "Herbicides", "Detergent Metabolites",
    "Antioxidants", "Antimicrobial Disinfectants",
    "Fire Retardants"
  ))

  plot_eps <- plot_tox_endpoints(chemical_summary,
    "Chemical",
    top_num = 5
  )

  expect_equal(
    tail(levels(plot_eps$data$endPoint), 5),
    c(
      "LTEA_HepaRG_EZR",
      "NVS_NR_mERa",
      "LTEA_HepaRG_HSPA1A",
      "LTEA_HepaRG_GSTM3",
      "CLD_HMGCS2_48hr"
    )
  )

  plot_stack <- plot_tox_stacks(chemical_summary,
    category = "Chemical",
    chem_site = tox_list$chem_site,
    top_num = 5
  )

  expect_equal(
    levels(plot_stack$data$category),
    c(
      "Bisphenol A", "4-Nonylphenol, branched",
      "Atrazine",
      "Pentachlorophenol", "Metolachlor",
      "Others (42)"
    )
  )

  plot_heat <- plot_tox_heatmap(chemical_summary,
    category = "Chemical",
    chem_site = tox_list$chem_site
  )

  expect_equal(
    tail(levels(plot_heat$data$chnm), 5),
    c(
      "Bromacil", 
      "Metalaxyl",
      "Metolachlor",
      "Atrazine",
      "Pentachlorophenol"
    )
  )
})
