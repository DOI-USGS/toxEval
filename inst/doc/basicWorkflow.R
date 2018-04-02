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


## ----boxplots1, warning=FALSE--------------------------------------------
plot_tox_boxplots(chemicalSummary, "Biological")   
plot_tox_boxplots(chemicalSummary, "Chemical Class")
plot_tox_boxplots(chemicalSummary, "Chemical") 

## ----stackplots1, warning=FALSE------------------------------------------
plot_tox_stacks(chemicalSummary, tox_list$chem_site, "Biological")   
plot_tox_stacks(chemicalSummary, tox_list$chem_site, "Chemical Class")
plot_tox_stacks(chemicalSummary, tox_list$chem_site, "Chemical", include_legend = FALSE) 

## ----heat, warning=FALSE-------------------------------------------------
plot_tox_heatmap(chemicalSummary, 
                  tox_list$chem_site, 
                  category = "Biological",
                  manual_remove = "Undefined")
plot_tox_heatmap(chemicalSummary, tox_list$chem_site, category = "Chemical Class")
plot_tox_heatmap(chemicalSummary, tox_list$chem_site, category = "Chemical")

## ----endpoints, warning=FALSE--------------------------------------------
plot_tox_endpoints(chemicalSummary, filterBy = "Cell Cycle")
plot_tox_endpoints(chemicalSummary, category = "Chemical Class", filterBy = "PAHs")
plot_tox_endpoints(chemicalSummary, category = "Chemical", filterBy = "Atrazine")

## ----table_tox_rank, warning=FALSE---------------------------------------
library(DT)
options(DT.options = list(pageLength = 5))

table_tox_rank(chemicalSummary, category = "Biological")
table_tox_rank(chemicalSummary, category = "Chemical Class")
table_tox_rank(chemicalSummary, category = "Chemical")

## ----table_tox_sum, warning=FALSE----------------------------------------
table_tox_sum(chemicalSummary, category = "Biological")
table_tox_sum(chemicalSummary, category = "Chemical Class")
table_tox_sum(chemicalSummary, category = "Chemical")

