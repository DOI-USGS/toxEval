## ----setup, include=FALSE---------------------------------
library(knitr)
library(rmarkdown)
options(continue=" ")
options(width=60)
knitr::opts_chunk$set(echo = TRUE,
                      warning = FALSE,
                      message = FALSE,
                      fig.height = 7,
                      fig.width = 7)

## ----getChems---------------------------------------------
library(toxEval)
library(dplyr)

path_to_tox <-  system.file("extdata", package="toxEval")
file_name <- "OWC_data_fromSup.xlsx"
full_path <- file.path(path_to_tox, file_name)

tox_list <- create_toxEval(full_path)

ACClong <- get_ACC(tox_list$chem_info$CAS)
ACClong <- remove_flags(ACClong)

cleaned_ep <- clean_endPoint_info(endPointInfo)
filtered_ep <- filter_groups(cleaned_ep)

chemicalSummary <- get_chemical_summary(tox_list, ACClong, filtered_ep)



## ----getAOPs----------------------------------------------
AOP_crosswalk <- read.csv(file.path(path_to_tox, "AOP_crosswalk.csv"), stringsAsFactors = FALSE)

## ----sumAOP-----------------------------------------------
AOP_summaries <- chemicalSummary %>%
  left_join(select(AOP_crosswalk,
                   endPoint=Component.Endpoint.Name,
                   AOP_id = AOP..,
                   AOP_title = AOP.Title), by="endPoint") %>%
  filter(!is.na(AOP_id)) %>%
  group_by(shortName,date,AOP_id, AOP_title) %>%
  summarise(sumEAR = sum(EAR, na.rm = TRUE)) %>%
  group_by(shortName, AOP_id, AOP_title) %>%
  summarise(maxEAR = max(sumEAR, na.rm = TRUE)) %>%
  group_by(AOP_id, AOP_title) %>%
  summarise(sum_of_maxes = sum(maxEAR)) %>%
  arrange(desc(sum_of_maxes)) 

threshold <- 1


kable(filter(AOP_summaries, sum_of_maxes> threshold))


## ----sumKE------------------------------------------------

KE_summaries <- chemicalSummary %>%
  left_join(select(AOP_crosswalk,
                   endPoint=Component.Endpoint.Name,
                   KE_id = KE.,
                   KE_title = Key.Event.Name), by="endPoint") %>%
  filter(!is.na(KE_id)) %>%
  group_by(shortName,date,KE_id, KE_title) %>%
  summarise(sumEAR = sum(EAR, na.rm = TRUE)) %>%
  group_by(shortName, KE_id, KE_title) %>%
  summarise(maxEAR = max(sumEAR, na.rm = TRUE)) %>%
  group_by(KE_id, KE_title) %>%
  summarise(sum_of_maxes = sum(maxEAR)) %>%
  arrange(desc(sum_of_maxes)) 

threshold <- 1

kable(filter(KE_summaries, sum_of_maxes> threshold))


## ----otherStuff-------------------------------------------

non_AOP_summaries <- chemicalSummary %>%
  anti_join(select(AOP_crosswalk,
                   endPoint=Component.Endpoint.Name), by="endPoint") %>%
  group_by(shortName,date, endPoint) %>%
  summarise(sumEAR = sum(EAR, na.rm = TRUE)) %>%
  group_by(shortName, endPoint) %>%
  summarise(maxEAR = max(sumEAR, na.rm = TRUE)) %>%
  group_by(endPoint) %>%
  summarise(sum_of_maxes = sum(maxEAR)) %>%
  arrange(desc(sum_of_maxes)) 

threshold <- 1

kable(filter(non_AOP_summaries, sum_of_maxes> threshold))



