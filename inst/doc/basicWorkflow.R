## ----startup, message=FALSE----------------------------------------------
library(toxEval)
path_to_tox <-  system.file("extdata", package="toxEval")
file_name <- "OWC_data_fromSup.xlsx"
full_path <- file.path(path_to_tox, file_name)

tox_list <- create_toxEval(full_path)


## ----chemical_summary----------------------------------------------------
ACC <- get_ACC(tox_list$chem_info$CAS)
ACC <- remove_flags(ACC)

cleaned_ep <- clean_endPoint_info(end_point_info)
filtered_ep <- filter_groups(cleaned_ep)

chemical_summary <- get_chemical_summary(tox_list, ACC, filtered_ep)

## ----eval=FALSE----------------------------------------------------------
#  names(end_point_info)

## ------------------------------------------------------------------------
cleaned_ep <- clean_endPoint_info(end_point_info)

filtered_ep <- filter_groups(cleaned_ep,
              groupCol = "intended_target_family",
              assays = c("ATG","NVS", "OT", "TOX21", 
                         "CEETOX", "APR", "CLD", "TANGUAY",
                         "NHEERL_PADILLA","NCCT_SIMMONS", "ACEA"),
              remove_groups = c("Background Measurement",
                                "Undefined"))

## ------------------------------------------------------------------------
unique(cleaned_ep$intended_target_family)

## ----eval=FALSE----------------------------------------------------------
#  unique(end_point_info$intended_target_family_sub)

## ----boxplots1, warning=FALSE, message=FALSE-----------------------------
plot_tox_boxplots(chemical_summary, "Biological")

# Other options:
# plot_tox_boxplots(chemical_summary, "Chemical Class")
# plot_tox_boxplots(chemical_summary, "Chemical") 

## ----plot_box_thres, warning=FALSE, message=FALSE------------------------
plot_tox_boxplots(chemical_summary, 
                  category = "Biological",
                  hit_threshold = 0.001)


## ----filtersiteBox, message=FALSE, warning=FALSE-------------------------
library(dplyr)

maumee <- filter(chemical_summary, shortName == "Maumee")
maumee_site <- filter(tox_list$chem_site, `Short Name` == "Maumee")

## ----maumeePlot, message=FALSE, warning=FALSE----------------------------
library(ggplot2)

plot_tox_boxplots(maumee, "Biological",title = maumee_site$Fullname[1])


## ----stackplots1, warning=FALSE, fig.width=10----------------------------
plot_tox_stacks(chemical_summary, 
                chem_site = tox_list$chem_site, 
                category =  "Biological")

# More options:
# plot_tox_stacks(chemical_summary, 
#                 chem_site = tox_list$chem_site, 
#                 category = "Chemical Class")
# plot_tox_stacks(chemical_summary, 
#                 chem_site = tox_list$chem_site, 
#                 category = "Chemical", include_legend = FALSE)

## ----siteStacks, message=FALSE, warning=FALSE, fig.width=10--------------
plot_tox_stacks(maumee, maumee_site,"Biological", title = maumee_site$Fullname[1])


## ----heat, warning=FALSE, fig.width=10-----------------------------------
plot_tox_heatmap(chemical_summary, 
                 chem_site = tox_list$chem_site, 
                 category = "Biological")
# More options:
# plot_tox_heatmap(chemical_summary, 
#                  chem_site = tox_list$chem_site, 
#                  category = "Chemical Class")
# plot_tox_heatmap(chemical_summary, 
#                  chem_site = tox_list$chem_site, 
#                  category = "Chemical")

## ----endpoints, warning=FALSE--------------------------------------------
plot_tox_endpoints(chemical_summary, 
                    category = "Biological", 
                    filterBy = "Cell Cycle")

# More options:
# plot_tox_endpoints(chemical_summary,   
#                    category = "Chemical Class", 
#                    filterBy = "PAHs")
# plot_tox_endpoints(chemical_summary, 
#                    category = "Chemical", 
#                    filterBy = "Atrazine")

## ----boxDetails, echo=FALSE----------------------------------------------
sample_df <- data.frame(parameter = "test",
                        values = rnorm(510, mean = 500, sd = 300))

#Make sure there are some outliers:
sample_df$values[1] <- 1000
sample_df$values[2] <- 900
sample_df$values[3] <- -200
sample_df$values[4] <- -600

quartiles <- as.numeric(quantile(sample_df$values))
IQR <- diff(quartiles[c(2,4)])

upper_whisker <- max(sample_df$values[sample_df$values <
                           (quartiles[4] + 1.5 * IQR)])
lower_whisker <- min(sample_df$values[sample_df$values >
                           (quartiles[2] - 1.5 * IQR)])
  
upper_dots <- sample_df$values[sample_df$values > 
                            quartiles[4] + 1.5*IQR]
lower_dots <- sample_df$values[sample_df$values < 
                            quartiles[2] - 1.5*IQR]

label_df <- data.frame(parameter = "explaination",
                       values = c(quartiles,
                                  upper_whisker,
                                  lower_whisker),
                       names = c("Minimum", 
                                 "Q1 = 25th percentile", "Median",
                                 "Q3 = 75th percentile", "Maximum","",""))


ggplot() +
  geom_hline(yintercept = label_df$values,linetype="dashed", color = "lightgrey") +
  geom_boxplot(data = sample_df, aes(x = parameter, y=values), width = 0.2, fill = "steelblue") +  theme_bw() +
  geom_text(data = label_df, aes(x = 0.75, y = values, label = names), vjust = 0.4) +
  geom_segment(aes(x = 1.15, xend = 1.15, 
                   y = label_df$values[2], yend = label_df$values[4]),
               arrow = arrow(length = unit(0.2, "cm"), ends = "both")) +
  geom_text(aes(x = 1.25, y = label_df$values[3]), 
            label = "IQR = Q3-Q1", hjust = 0.2, vjust = 0.4) +
  geom_text(aes(x = c(1.1,1.1), y = label_df$values[6:7]), 
            label = c(bquote(italic(max)*"(y[y < (Q3 + 1.5*IQR)])"),
                      bquote(italic(min)*"(y[y > (Q1 - 1.5*IQR)])")), 
            hjust = 0, vjust = 0.4, parse = TRUE) +
  geom_text(aes(x = c(1.1), y = mean(upper_dots)), 
            label = "upper dots > Q3 + 1.5*IQR", hjust = 0, vjust = 0.4) +
  geom_text(aes(x = c(1.1), y = mean(lower_dots)), 
            label = "lower dots < Q1 - 1.5*IQR", hjust = 0, vjust = 0.4) +
  ylab("") + xlab("") + 
  labs(title = "ggplot2 boxplot explanation", 
       caption = bquote(italic("Q = quartiles, IQR = Interquartile range, y = data"))) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())


## ----ggsave1, eval=FALSE-------------------------------------------------
#  
#  ep_plot <- plot_tox_endpoints(chemical_summary,
#                                category = "Biological",
#                                filterBy = "Cell Cycle")
#  
#  # To save a png:
#  ggsave(ep_plot, file = "ep_plot.png")
#  
#  # To save a pdf:
#  ggsave(ep_plot, file = "ep_plot.pdf")

## ----rank_sites_DT, warning=FALSE----------------------------------------
library(DT)
options(DT.options = list(pageLength = 5))

rank_df <- rank_sites(chemical_summary, 
                      category = "Biological",
                      hit_threshold = 0.1)

rank_sites_DT(chemical_summary, 
              category = "Biological",
              hit_threshold = 0.1)

## ----rank_sites_DT_site, warning=FALSE-----------------------------------
rank_sites_DT(maumee, category = "Biological")

## ----hits_summary_DT, warning=FALSE--------------------------------------

hit_df <- hits_summary(chemical_summary,
                       category = "Biological",
                       hit_threshold = 0.1 )

hits_summary_DT(chemical_summary, 
                category = "Biological",
                hit_threshold = 0.1)

## ----hits_summary_DT_site, warning=FALSE---------------------------------
hits_summary_DT(maumee, category = "Biological")

## ----endpoint_hits_DT, warning=FALSE-------------------------------------

ep_hits <- endpoint_hits(chemical_summary, 
                         category = "Biological", 
                         hit_threshold = 0.1)

endpoint_hits_DT(chemical_summary, 
                 category = "Biological",
                 hit_threshold = 0.1)


## ----endpoint_hits_DT_site, warning=FALSE--------------------------------
endpoint_hits_DT(maumee, category = "Biological")

## ----hits_by_groupings_DT, warning=FALSE---------------------------------
site_df <- hits_by_groupings(chemical_summary, 
                             category = "Chemical Class",
                             hit_threshold = 0.1)

hits_by_groupings_DT(chemical_summary, 
                     category = "Chemical Class",
                     hit_threshold = 0.1)

## ----hits_by_groupings_DT_site, warning=FALSE----------------------------
hits_by_groupings_DT(maumee, category = "Chemical Class")

## ----makeMap, warning=FALSE, message=FALSE-------------------------------
make_tox_map(chemical_summary, 
             chem_site = tox_list$chem_site, 
             category = "Biological")
# More options:
# make_tox_map(chemical_summary, 
#              chem_site = tox_list$chem_site, 
#              category = "Chemical Class")
# make_tox_map(chemical_summary, 
#              chem_site = tox_list$chem_site, 
#              category = "Chemical")


## ----clean---------------------------------------------------------------
#Trim some names:
levels(chemical_summary$Class)[levels(chemical_summary$Class) == "Antimicrobial Disinfectants"] <- "Antimicrobial"
levels(chemical_summary$Class)[levels(chemical_summary$Class) == "Detergent Metabolites"] <- "Detergent"
levels(chemical_summary$Class)[levels(chemical_summary$Class) == "Flavors and Fragrances"] <- "Flavor/Fragrance"

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


## ----fig.width=10--------------------------------------------------------
summary_with_levels <- get_chemical_summary(tox_list,
                                            ACC,
                                            filtered_ep)

plot_tox_stacks(summary_with_levels, tox_list$chem_site, "Biological")
plot_tox_heatmap(summary_with_levels, tox_list$chem_site, "Biological")

