## ----startup, message=FALSE----------------------------------------------
library(toxEval)
path_to_tox <-  system.file("extdata", package="toxEval")
file_name <- "OWC_data_fromSup.xlsx"
full_path <- file.path(path_to_tox, file_name)

tox_list <- create_toxEval(full_path)


## ----chemicalSummary-----------------------------------------------------
ACClong <- get_ACC(tox_list$chem_info$CAS)
ACClong <- remove_flags(ACClong)

cleaned_ep <- clean_endPoint_info(endPointInfo)
filtered_ep <- filter_groups(cleaned_ep)

chemicalSummary <- get_chemical_summary(tox_list, ACClong, filtered_ep)

## ----eval=FALSE----------------------------------------------------------
#  names(endPointInfo)

## ----clean---------------------------------------------------------------
#Trim some names:
levels(chemicalSummary$Class)[levels(chemicalSummary$Class) == "Antimicrobial Disinfectants"] <- "Antimicrobial"
levels(chemicalSummary$Class)[levels(chemicalSummary$Class) == "Detergent Metabolites"] <- "Detergent"
levels(chemicalSummary$Class)[levels(chemicalSummary$Class) == "Flavors and Fragrances"] <- "Flavor/Fragrance"

## ----sites---------------------------------------------------------------
#Ordering the sites to flow "downstream" of the Great Lakes:
sitesOrdered <- c("StLouis","Nemadji","WhiteWI","Bad","Montreal",
  "PresqueIsle","Ontonagon","Sturgeon","Tahquamenon","Burns",
  "IndianaHC","StJoseph","PawPaw","Kalamazoo","GrandMI",
  "Milwaukee","Muskegon","WhiteMI","PereMarquette","Manitowoc",
  "Manistee","Fox","Oconto","Peshtigo","Menominee",
  "Indian","Cheboygan","Ford","Escanaba","Manistique",
  "ThunderBay","AuSable","Rifle","Saginaw","BlackMI",
  "Clinton","Rouge","HuronMI","Raisin","Maumee",
  "Portage","Sandusky","HuronOH","Vermilion","BlackOH",
  "Rocky","Cuyahoga","GrandOH","Cattaraugus","Tonawanda",
  "Genesee","Oswego","BlackNY","Oswegatchie","Grass",
  "Raquette","StRegis")

 tox_list$chem_site$`Short Name` <- factor(tox_list$chem_site$`Short Name`,
               levels = sitesOrdered)

lakes_ordered <- c("Lake Superior",
                  "Lake Michigan",
                  "Lake Huron",
                  "Lake Erie",
                  "Lake Ontario")

tox_list$chem_site$site_grouping <- factor(tox_list$chem_site$site_grouping,
               levels=lakes_ordered)


## ----boxplots1, warning=FALSE, message=FALSE-----------------------------

bio_box <- plot_tox_boxplots(chemicalSummary, "Biological")

# The graph can be plotted without these additional lines,
# but they allow the labels to look nicer:
gb <- ggplot2::ggplot_build(bio_box)
gt <- ggplot2::ggplot_gtable(gb)
gt$layout$clip[gt$layout$name=="panel"] <- "off"
grid::grid.draw(gt)
# Other options:
# plot_tox_boxplots(chemicalSummary, "Chemical Class")
# plot_tox_boxplots(chemicalSummary, "Chemical") 

## ----filtersiteBox, message=FALSE, warning=FALSE-------------------------
library(dplyr)

maumee <- filter(chemicalSummary, shortName == "Maumee")
maumee_site <- filter(tox_list$chem_site, `Short Name` == "Maumee")

## ----maumeePlot, message=FALSE, warning=FALSE----------------------------
library(ggplot2)

maumee_plot <- plot_tox_boxplots(maumee, "Biological",title = maumee_site$Fullname[1])

gb <- ggplot2::ggplot_build(maumee_plot)
gt <- ggplot2::ggplot_gtable(gb)
gt$layout$clip[gt$layout$name=="panel"] <- "off"
grid::grid.draw(gt)


## ----stackplots1, warning=FALSE, fig.width=10--------------
stack_plot <- plot_tox_stacks(chemicalSummary, tox_list$chem_site, "Biological")

gb <- ggplot2::ggplot_build(stack_plot)
gt <- ggplot2::ggplot_gtable(gb)
gt$layout$clip[gt$layout$name=="panel-1-1"] = "off"
grid::grid.draw(gt)

# More options:
# plot_tox_stacks(chemicalSummary, tox_list$chem_site, "Chemical Class")
# plot_tox_stacks(chemicalSummary, tox_list$chem_site, "Chemical", include_legend = FALSE) 

## ----siteStacks, message=FALSE, warning=FALSE, fig.width=10----
maumee_plot_stack <- plot_tox_stacks(maumee, maumee_site,"Biological", title = maumee_site$Fullname[1])

gb <- ggplot2::ggplot_build(maumee_plot_stack)
gt <- ggplot2::ggplot_gtable(gb)
gt$layout$clip[gt$layout$name=="panel"] <- "off"
grid::grid.draw(gt)


## ----heat, warning=FALSE, fig.width=10---------------------
plot_tox_heatmap(chemicalSummary, 
                  tox_list$chem_site, 
                  category = "Biological")
# More options:
# plot_tox_heatmap(chemicalSummary, tox_list$chem_site, category = "Chemical Class")
# plot_tox_heatmap(chemicalSummary, tox_list$chem_site, category = "Chemical")

## ----endpoints, warning=FALSE------------------------------
ep_plot <- plot_tox_endpoints(chemicalSummary, filterBy = "Cell Cycle")

gb <- ggplot2::ggplot_build(ep_plot)
gt <- ggplot2::ggplot_gtable(gb)
gt$layout$clip[gt$layout$name=="panel"] <- "off"
grid::grid.draw(gt)

# More options:
# plot_tox_endpoints(chemicalSummary, category = "Chemical Class", filterBy = "PAHs")
# plot_tox_endpoints(chemicalSummary, category = "Chemical", filterBy = "Atrazine")

## ----table_tox_rank, warning=FALSE-------------------------
library(DT)
options(DT.options = list(pageLength = 5))

table_tox_rank(chemicalSummary, category = "Biological")
# More options:
# table_tox_rank(chemicalSummary, category = "Chemical Class")
# table_tox_rank(chemicalSummary, category = "Chemical")

## ----table_tox_rank_site, warning=FALSE--------------------
table_tox_rank(maumee, category = "Biological")

## ----table_tox_sum, warning=FALSE--------------------------
table_tox_sum(chemicalSummary, category = "Biological")
# More options:
# table_tox_sum(chemicalSummary, category = "Chemical Class")
# table_tox_sum(chemicalSummary, category = "Chemical")

## ----table_tox_sum_site, warning=FALSE---------------------
table_tox_sum(maumee, category = "Biological")

## ----table_endpoint_hits, warning=FALSE--------------------
table_endpoint_hits(chemicalSummary, category = "Biological")
# More options:
# table_endpoint_hits(chemicalSummary, category = "Chemical Class")
# table_endpoint_hits(chemicalSummary, category = "Chemical")

## ----table_endpoint_hits_site, warning=FALSE---------------
table_endpoint_hits(maumee, category = "Biological")

## ----table_tox_endpoint, warning=FALSE---------------------
table_tox_endpoint(chemicalSummary, category = "Chemical Class")
# More options:
# table_tox_endpoint(chemicalSummary, category = "Biological")
# table_tox_endpoint(chemicalSummary, category = "Chemical")

## ----table_tox_endpoint_site, warning=FALSE----------------
table_tox_endpoint(maumee, category = "Chemical Class")

## ----makeMap, warning=FALSE, message=FALSE-----------------
make_tox_map(chemicalSummary, tox_list$chem_site, "Biological")
# More options:
# make_tox_map(chemicalSummary, tox_list$chem_site, "Chemical Class")
# make_tox_map(chemicalSummary, tox_list$chem_site, "Chemical") 


