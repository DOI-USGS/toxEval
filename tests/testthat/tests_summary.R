context("Chemical summary tests")

library(readxl)
path_to_tox <-  system.file("extdata", package="toxEval")
file_name <- "OWC_data_fromSup.xlsx"
full_path <- file.path(path_to_tox, file_name)

chem_data <- read_excel(full_path, sheet = "Data")
chem_info <- read_excel(full_path, sheet = "Chemicals")
chem_site <- read_excel(full_path, sheet = "Sites")
ACClong <- get_ACC(chem_info$CAS)
ACClong <- remove_flags(ACClong)

cleaned_ep <- clean_endPoint_info(endPointInfo)
filtered_ep <- filter_groups(cleaned_ep)

chemicalSummary <- get_chemical_summary(ACClong,
                                        filtered_ep,
                                        chem_data,
                                        chem_site,
                                        chem_info)

test_that("Calculating summaries", {
  testthat::skip_on_cran()
 
  expect_true(all(names(chemicalSummary) %in% c("casrn","chnm","endPoint","site","date",
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
  expect_equal(bioPlot$layers[[1]]$geom_params$outlier.shape, 19)
  expect_equal(bioPlot$layers[[1]]$aes_params$fill, "steelblue")
  
  classPlot <- plot_tox_boxplots(chemicalSummary, 
                               category = "Chemical Class")

  expect_true(all(names(classPlot$data) %in% c("site","category","meanEAR")))
  expect_equal(classPlot$layers[[1]]$geom_params$outlier.shape, 19)
  expect_equal(classPlot$layers[[1]]$aes_params$fill, "steelblue")
  
  chemPlot <- suppressWarnings(plot_tox_boxplots(chemicalSummary, 
                                 category = "Chemical"))
  
  expect_true(all(names(chemPlot$data) %in% c("site","chnm","Class","maxEAR")))
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
  expect_true(all(names(chemHeatPlot$data) %in% c("site","chnm","Class","maxEAR",
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
                                   category = "Chemical")
  expect_true(all(names(chemStackPlot$data) %in% c("site","category","meanEAR",
                                                  "site_grouping","Short Name")))
  
})


test_that("Internal table functions", {
  testthat::skip_on_cran()
  
  statStuff <- statsOfColumns(chemicalSummary, "Biological", hit_threshold = 0.1, mean_logic = FALSE)
  expect_equal(nrow(statStuff), nrow(chem_site))
  
  expect_true(all(c("Cell Cycle freq", "Cell Cycle maxEAR") %in% names(statStuff)))
  
  groupStuff <- statsOfGroup(chemicalSummary, "Biological", hit_threshold = 1)
  
  expect_true(all(unique(groupStuff$site) %in% chem_site$`Short Name`))
})
