context("Utilities")

test_that("Labels", {
  testthat::skip_on_cran()

  lab1 <- toxEval:::fancyLabels("Chemical Class",
    mean_logic = FALSE, sum_logic = TRUE, sep = TRUE,
    single_site = FALSE, include_site = TRUE
  )

  expect_equal(lab1[["caption"]], bquote(italic("i = chemicals in a specified class, j = samples, k = sites")))
  expect_equal(lab1[["y_label"]], bquote(italic("max") ~ group("[", group("(", sum(" " * EAR["[" * i * "]"]), ")")["[" * j * "]"], "]")["[" * k * "]"]))

  lab2 <- toxEval:::fancyLabels("Chemical",
    mean_logic = TRUE, sum_logic = TRUE, sep = TRUE,
    single_site = FALSE, include_site = TRUE
  )

  expect_equal(lab2[["caption"]], bquote(italic("i = chemicals, j = samples, k = sites")))
  expect_equal(lab2[["y_label"]], bquote(italic("mean") ~ group("[", group("(", sum(" " *
    EAR["[" * i * "]"]), ")")["[" * j * "]"], "]")["[" * k * "]"]))
  lab3 <- toxEval:::fancyLabels("Biological",
    mean_logic = TRUE, sum_logic = FALSE, sep = TRUE,
    single_site = FALSE, include_site = TRUE
  )

  expect_equal(lab3[["caption"]], bquote(italic("i = chemicals in a specified grouping, j = samples, k = sites")))
  expect_equal(lab3[["y_label"]], bquote(italic("mean") * group("[", italic(max) * group("(", EAR["[" * i * "]"], ")")["[" * j * "]"], "]")["[" * k * "]"]))

  nums1 <- toxEval:::fancyNumbers(c(10, 100, 1000, 10000))
  expect_equal(nums1[[1]], 10)
  expect_equal(as.character(nums1[[4]]), c("^", "10", "4"))

  nums2 <- toxEval:::fancyNumbers2(c(10, 100, 1000, 10000))
  expect_equal(nums2, c("10", "100", "1000", "10000"))

  nums3 <- toxEval:::prettyLogs(c(0.01, 0.5, 1, 20, 300))
  expect_equal(nums3, c(0.01, 0.1, 1, 10, 100, 1000))
})


