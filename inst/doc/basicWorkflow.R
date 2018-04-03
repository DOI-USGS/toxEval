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
library(grid)
plot_tox_boxplots(chemicalSummary, "Biological")   
grid.text("# Sites:", 
          x = unit(.22, "npc"), 
          y = unit(.995, "npc"), gp=gpar(fontsize=7))
# Other options:
# plot_tox_boxplots(chemicalSummary, "Chemical Class")
# plot_tox_boxplots(chemicalSummary, "Chemical") 

## ----siteBox, message=FALSE, warning=FALSE-------------------------------
library(dplyr)
maumee <- filter(chemicalSummary, shortName == "Maumee")

plot_tox_boxplots(maumee, "Biological")
grid.text("# EndPoints:", 
          x = unit(.22, "npc"), 
          y = unit(.995, "npc"), gp=gpar(fontsize=7))

## ----stackplots1, warning=FALSE, fig.width=10----------------------------
plot_tox_stacks(chemicalSummary, tox_list$chem_site, "Biological")
grid.text("# Detections:", 
          x = unit(.05, "npc"), 
          y = unit(.205, "npc"), gp=gpar(fontsize=7))
# More options:
# plot_tox_stacks(chemicalSummary, tox_list$chem_site, "Chemical Class")
# plot_tox_stacks(chemicalSummary, tox_list$chem_site, "Chemical", include_legend = FALSE) 

## ----siteStacks, message=FALSE, warning=FALSE, fig.width=10--------------

maumee_site <- filter(tox_list$chem_site, `Short Name` == "Maumee")
plot_tox_stacks(maumee, maumee_site,"Biological")


## ----heat, warning=FALSE, fig.width=10-----------------------------------
plot_tox_heatmap(chemicalSummary, 
                  tox_list$chem_site, 
                  category = "Biological")
# More options:
# plot_tox_heatmap(chemicalSummary, tox_list$chem_site, category = "Chemical Class")
# plot_tox_heatmap(chemicalSummary, tox_list$chem_site, category = "Chemical")

## ----endpoints, warning=FALSE--------------------------------------------
plot_tox_endpoints(chemicalSummary, filterBy = "Cell Cycle")
grid.text("# Detections:", 
          x = unit(.38, "npc"), 
          y = unit(.995, "npc"), gp=gpar(fontsize=7))
# More options:
# plot_tox_endpoints(chemicalSummary, category = "Chemical Class", filterBy = "PAHs")
# plot_tox_endpoints(chemicalSummary, category = "Chemical", filterBy = "Atrazine")

## ----table_tox_rank, warning=FALSE---------------------------------------
library(DT)
options(DT.options = list(pageLength = 5))

table_tox_rank(chemicalSummary, category = "Biological")
# More options:
# table_tox_rank(chemicalSummary, category = "Chemical Class")
# table_tox_rank(chemicalSummary, category = "Chemical")

## ----table_tox_rank_site, warning=FALSE----------------------------------
table_tox_rank(maumee, category = "Biological")

## ----table_tox_sum, warning=FALSE----------------------------------------
table_tox_sum(chemicalSummary, category = "Biological")
# More options:
# table_tox_sum(chemicalSummary, category = "Chemical Class")
# table_tox_sum(chemicalSummary, category = "Chemical")

## ----table_tox_sum_site, warning=FALSE-----------------------------------
table_tox_sum(maumee, category = "Biological")

## ----table_endpoint_hits, warning=FALSE----------------------------------
table_endpoint_hits(chemicalSummary, category = "Biological")
# More options:
# table_endpoint_hits(chemicalSummary, category = "Chemical Class")
# table_endpoint_hits(chemicalSummary, category = "Chemical")

## ----table_endpoint_hits_site, warning=FALSE-----------------------------
table_endpoint_hits(maumee, category = "Biological")

## ----table_tox_endpoint, warning=FALSE-----------------------------------
table_tox_endpoint(chemicalSummary, category = "Biological")
# More options:
# table_tox_endpoint(chemicalSummary, category = "Chemical Class")
# table_tox_endpoint(chemicalSummary, category = "Chemical")

## ----table_tox_endpoint_site, warning=FALSE------------------------------
table_tox_endpoint(maumee, category = "Biological")

