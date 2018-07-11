path_to_excel <- "neonic_bench.xlsx"
tox_list <- create_toxEval(path_to_excel)
chemicalSummary <- get_chemical_summary(tox_list)

plot_tox_boxplots(chemicalSummary, category = "Chemical")
plot_tox_boxplots(chemicalSummary, category = "Chemical Class")

plot_tox_boxplots(chemicalSummary, category = "Biological")

