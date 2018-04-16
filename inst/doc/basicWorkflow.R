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

## ----plot_box_thres, warning=FALSE, message=FALSE------------------------
bio_box_thresh <- plot_tox_boxplots(chemicalSummary, 
                             category = "Biological",
                             hit_threshold = 0.001)

# The graph can be plotted without these additional lines,
# but they allow the labels to look nicer:
gb <- ggplot2::ggplot_build(bio_box_thresh)
gt <- ggplot2::ggplot_gtable(gb)
gt$layout$clip[gt$layout$name=="panel"] <- "off"
grid::grid.draw(gt)


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


## ----stackplots1, warning=FALSE, fig.width=10----------------------------
stack_plot <- plot_tox_stacks(chemicalSummary, 
                              chem_site = tox_list$chem_site, 
                              category =  "Biological")

gb <- ggplot2::ggplot_build(stack_plot)
gt <- ggplot2::ggplot_gtable(gb)
gt$layout$clip[gt$layout$name=="panel-1-1"] = "off"
grid::grid.draw(gt)

# More options:
# plot_tox_stacks(chemicalSummary, 
#                 chem_site = tox_list$chem_site, 
#                 category = "Chemical Class")
# plot_tox_stacks(chemicalSummary, 
#                 chem_site = tox_list$chem_site, 
#                 category = "Chemical", include_legend = FALSE)

## ----siteStacks, message=FALSE, warning=FALSE, fig.width=10--------------
maumee_plot_stack <- plot_tox_stacks(maumee, maumee_site,"Biological", title = maumee_site$Fullname[1])

gb <- ggplot2::ggplot_build(maumee_plot_stack)
gt <- ggplot2::ggplot_gtable(gb)
gt$layout$clip[gt$layout$name=="panel"] <- "off"
grid::grid.draw(gt)


## ----heat, warning=FALSE, fig.width=10-----------------------------------
plot_tox_heatmap(chemicalSummary, 
                 chem_site = tox_list$chem_site, 
                 category = "Biological")
# More options:
# plot_tox_heatmap(chemicalSummary, 
#                  chem_site = tox_list$chem_site, 
#                  category = "Chemical Class")
# plot_tox_heatmap(chemicalSummary, 
#                  chem_site = tox_list$chem_site, 
#                  category = "Chemical")

## ----endpoints, warning=FALSE--------------------------------------------
ep_plot <- plot_tox_endpoints(chemicalSummary, 
                              category = "Biological", 
                              filterBy = "Cell Cycle")

gb <- ggplot2::ggplot_build(ep_plot)
gt <- ggplot2::ggplot_gtable(gb)
gt$layout$clip[gt$layout$name=="panel"] <- "off"
grid::grid.draw(gt)

# More options:
# plot_tox_endpoints(chemicalSummary,   
#                    category = "Chemical Class", 
#                    filterBy = "PAHs")
# plot_tox_endpoints(chemicalSummary, 
#                    category = "Chemical", 
#                    filterBy = "Atrazine")

## ----clipping, eval=FALSE------------------------------------------------
#  gb <- ggplot2::ggplot_build(ep_plot)
#  gt <- ggplot2::ggplot_gtable(gb)
#  gt$layout$clip[gt$layout$name=="panel"] <- "off"
#  grid::grid.draw(gt)

## ----ggsave1, eval=FALSE-------------------------------------------------
#  
#  ep_plot <- plot_tox_endpoints(chemicalSummary,
#                                category = "Biological",
#                                filterBy = "Cell Cycle")
#  
#  # To save a png:
#  ggsave(ep_plot, file = "ep_plot.png")
#  
#  # To save a pdf:
#  ggsave(ep_plot, file = "ep_plot.pdf")

## ----basesave1, eval=FALSE-----------------------------------------------
#  
#  ep_plot <- plot_tox_endpoints(chemicalSummary,
#                                category = "Biological",
#                                filterBy = "Cell Cycle")
#  gb <- ggplot2::ggplot_build(ep_plot)
#  gt <- ggplot2::ggplot_gtable(gb)
#  gt$layout$clip[gt$layout$name=="panel"] <- "off"
#  
#  
#  # To save a png:
#  png("ep_plot.png", width = 1200, height = 1200, res = 142)
#  grid::grid.draw(gt)
#  dev.off()
#  
#  # To save a pdf:
#  pdf("ep_plot.pdf", width = 9, height = 11)
#  grid::grid.draw(gt)
#  dev.off()
#  
#  

## ----rank_sites_DT, warning=FALSE----------------------------------------
library(DT)
options(DT.options = list(pageLength = 5))

rank_df <- rank_sites(chemicalSummary, 
                      category = "Biological",
                      hit_threshold = 0.1)

rank_sites_DT(chemicalSummary, 
              category = "Biological",
              hit_threshold = 0.1)

## ----rank_sites_DT_site, warning=FALSE-----------------------------------
rank_sites_DT(maumee, category = "Biological")

## ----hits_summary_DT, warning=FALSE--------------------------------------

hit_df <- hits_summary(chemicalSummary,
                       category = "Biological",
                       hit_threshold = 0.1 )

hits_summary_DT(chemicalSummary, 
                category = "Biological",
                hit_threshold = 0.1)

## ----hits_summary_DT_site, warning=FALSE---------------------------------
hits_summary_DT(maumee, category = "Biological")

## ----endpoint_hits_DT, warning=FALSE-------------------------------------

ep_hits <- endpoint_hits(chemicalSummary, 
                         category = "Biological", 
                         hit_threshold = 0.1)

endpoint_hits_DT(chemicalSummary, 
                 category = "Biological",
                 hit_threshold = 0.1)


## ----endpoint_hits_DT_site, warning=FALSE--------------------------------
endpoint_hits_DT(maumee, category = "Biological")

## ----hits_by_groupings_DT, warning=FALSE---------------------------------
site_df <- hits_by_groupings(chemicalSummary, 
                             category = "Chemical Class",
                             hit_threshold = 0.1)

hits_by_groupings_DT(chemicalSummary, 
                     category = "Chemical Class",
                     hit_threshold = 0.1)

## ----hits_by_groupings_DT_site, warning=FALSE----------------------------
hits_by_groupings_DT(maumee, category = "Chemical Class")

## ----makeMap, warning=FALSE, message=FALSE-------------------------------
make_tox_map(chemicalSummary, 
             chem_site = tox_list$chem_site, 
             category = "Biological")
# More options:
# make_tox_map(chemicalSummary, 
#              chem_site = tox_list$chem_site, 
#              category = "Chemical Class")
# make_tox_map(chemicalSummary, 
#              chem_site = tox_list$chem_site, 
#              category = "Chemical")


