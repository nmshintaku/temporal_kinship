#######################################
#### Calculate Association Coefficients
#######################################

library(SocGen)
library(igraph)
library(dplyr)
library(tidyr)

indiv_covars <- read.csv("Outputs/individual_covariates.csv")
indiv_covars$entrydate <- as.Date(indiv_covars$entrydate)
indiv_covars$departdate <- as.Date(indiv_covars$departdate)

sightings <- read.csv("Outputs/sightings.csv")
pairwise <- read.csv("Outputs/pairwise_covariates.csv")

dates <- sort(unique(sightings$Observation.Date))

#Adding separate columns for Repro.ID and Age.ID
sightings <- sightings %>%
  mutate(Repro.ID = sub("^([^.]+\\.[^.]+).*", "\\1", Combined.ID)) %>%
  mutate(Age.ID = paste(sub("\\..*", "", Combined.ID), sub(".*\\.(\\w+)$", "\\1", Combined.ID), sep = "."))

#count repro sightings

rs_tab <- table(sightings$Combined.ID)
sightings$ReproSightings <- rs_tab

# select just genotyped individuals

relatedness <- read.csv("Raw Data/RelatednessEstimates_2024.csv")
genotypes <- unique(c(relatedness$ID1, relatedness$ID2))

indiv_covars$genotyped <- ifelse(indiv_covars$Dolphin.ID %in% genotypes, "Y", "N")

statuses <- sightings[,c("Dolphin.ID", "Combined.ID")] |> unique()

affil_females <- merge(statuses, indiv_covars, by = "Dolphin.ID", all.x = TRUE)

affil_females <- affil_females[which(affil_females$genotyped == "Y"), ]

#Adding separate columns for Repro.ID and Age.ID
affil_females <- affil_females %>%
  mutate(Repro.ID = sub("^([^.]+\\.[^.]+).*", "\\1", Combined.ID)) %>%
  mutate(Age.ID = paste(sub("\\..*", "", Combined.ID), sub(".*\\.(\\w+)$", "\\1", Combined.ID), sep = "."))

# mask the data so that association rates are only estimated in the timeframe where both members are alive
affil_mask <- schedulize(affil_females,
                         id = "Combined.ID",
                         start = "entrydate",
                         end = "departdate",
                         dates = dates,
                         format = "mask")

affil_sightings <- sightings[which(sightings$Combined.ID %in% affil_females$Combined.ID), ]

#Adding separate columns for Repro.ID and Age.ID
affil_sightings <- affil_sightings %>%
  mutate(Repro.ID = sub("^([^.]+\\.[^.]+).*", "\\1", Combined.ID)) %>%
  mutate(Age.ID = paste(sub("\\..*", "", Combined.ID), sub(".*\\.(\\w+)$", "\\1", Combined.ID), sep = "."))

#Calculating SRI
masked_network <- simple_ratio(sightings = affil_sightings,
                              group_variable = "Observation.ID",
                              dates = "Observation.Date",
                              IDs = "Combined.ID",
                              mask = affil_mask)

# For each female reprocat combination, sum up per reprocat and per relatedness cat 

# Create a network object

masked_network[is.nan(masked_network)] <- 0
masked_network[is.na(masked_network)] <- 0

network <- graph_from_adjacency_matrix(masked_network, mode = "undirected", weighted = TRUE)

# calculate strength (summed edge weights) per ID
tot_str <- strength(network)
tot_str_df <- data.frame(name = V(network)$name, strength = tot_str)
#write.csv(tot_str_df, "Outputs/total_strength.csv", row.names = FALSE)

#########
##Repro##
#########

# mask the data so that association rates are only estimated in the timeframe where both members are alive
affil_mask_repro <- schedulize(affil_females,
                         id = "Repro.ID",
                         start = "entrydate",
                         end = "departdate",
                         dates = dates,
                         format = "mask")

repro_mask_network <- simple_ratio(sightings = affil_sightings,
                               group_variable = "Observation.ID",
                               dates = "Observation.Date",
                               IDs = "Repro.ID",
                               mask = affil_mask_repro)

repro_mask_network[is.nan(repro_mask_network)] <- 0
repro_mask_network[is.na(repro_mask_network)] <- 0

repro_network <- graph_from_adjacency_matrix(repro_mask_network, mode = "undirected", weighted = TRUE)

str_repro <- strength(repro_network)
str_repro_df <- data.frame(name = V(repro_network)$name, strength = str_repro)
#write.csv(str_repro_df, "Outputs/repro_strength.csv", row.names = FALSE)

str_repro_df <- str_repro_df %>%
  separate(name, into = c("Dolphin.ID", "type"), sep = "\\.") %>%
  spread(key = type, value = strength) %>%
  replace_na(list(cyc = 0, lact = 0, preg = 0, juvenile = 0, unknown = 0)) %>%
  rename(lact.strength = lact, cyc.strength = cyc, preg.strength = preg, repro.juvenile.strength = juvenile, repro.unknown.strength = unknown)

#Not sure how to merge into affil_females
affil_females <- merge(affil_females, str_repro_df, by)

#CODE NOT WORKING 
merge_repro_str <- merge(sightings, str_repro_df, by.x = "Repro.ID", by.y = "name", all.x = TRUE)

#Setting vertex attributes to each reproductive category
repro_network <- set_vertex_attr(repro_network, "ReproCat",
                                 value = sightings$ReproCat[match(V(repro_network)$name, sightings$Repro.ID)])


#########
###Age###
#########

affil_mask_age <- schedulize(affil_females,
                               id = "Age.ID",
                               start = "entrydate",
                               end = "departdate",
                               dates = dates,
                               format = "mask")

age_mask_network <- simple_ratio(sightings = affil_sightings,
                                   group_variable = "Observation.ID",
                                   dates = "Observation.Date",
                                   IDs = "Age.ID",
                                   mask = affil_mask_age)

age_mask_network[is.nan(age_mask_network)] <- 0
age_mask_network[is.na(age_mask_network)] <- 0

age_network <- graph_from_adjacency_matrix(age_mask_network, mode = "undirected", weighted = TRUE)

str_age <- strength(age_network)
str_age_df <- data.frame(name = V(age_network)$name, strength = str_age)
#write.csv(str_age_df, "Outputs/age_strength.csv", row.names = FALSE)

str_age_df <- str_age_df %>%
  separate(name, into = c("Dolphin.ID", "type"), sep = "\\.") %>%
  spread(key = type, value = strength) %>%
  replace_na(list(juvenile = 0, adult = 0, elder = 0)) %>%
  rename(adult.strength = adult, juvenile.strength = juvenile, elder.strength = elder)

#Setting vertex attributes to each age class
age_network <- set_vertex_attr(age_network, "AgeClass",
                                 value = sightings$AgeClass[match(V(age_network)$name, sightings$Age.ID)])


#Combining all the strength calculations into one dataframe
all_strength <- merge(str_repro_df, str_age_df, by = "Dolphin.ID", all.x = TRUE)
write.csv(all_strength, "Outputs/all_strength.csv", row.names = FALSE)


