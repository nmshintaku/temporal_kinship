#02 Association Scratch script (Close Kin)

library(SocGen)
library(igraph)
library(dplyr)
library(tidyr)
library(ggplot2)

indiv_covars <- read.csv("Outputs/individual_covariates.csv")
indiv_covars$entrydate <- as.Date(indiv_covars$entrydate)
indiv_covars$departdate <- as.Date(indiv_covars$departdate)

sightings <- read.csv("Outputs/sightings.csv")
pairwise <- read.csv("Outputs/pairwise_covariates.csv")

#Removing NA in pairwise
pairwise <- pairwise %>%
  filter(!is.na(DyadML))

dates <- sort(unique(sightings$Observation.Date))

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

affil_sightings <- sightings[which(sightings$Combined.ID %in% affil_females$Combined.ID), ]

#Creating ReproCat column
extract_state <- function(Repro.ID) {
  return(sub(".*\\.", "", Repro.ID))
}

affil_females$ReproCat <- sapply(affil_females$Repro.ID, extract_state)

#Creating ReproCat column
extract_age <- function(Age.ID) {
  return(sub(".*\\.", "", Age.ID))
}

affil_females$AgeCat <- sapply(affil_females$Age.ID, extract_age)

# mask the data so that association rates are only estimated in the timeframe where both members are alive
affil_mask <- schedulize(affil_females,
                         id = "Combined.ID",
                         start = "entrydate",
                         end = "departdate",
                         dates = dates,
                         format = "mask")

#Calculating SRI
masked_network <- simple_ratio(sightings = affil_sightings,
                               group_variable = "Observation.ID",
                               dates = "Observation.Date",
                               IDs = "Combined.ID",
                               mask = affil_mask)

# Create a network object

masked_network[is.nan(masked_network)] <- 0
masked_network[is.na(masked_network)] <- 0

combined_list <- mat2dat(masked_network, value.name = "weight")

#Set kinship as the edge attribute 

#create list of all pairwise combinations of combined IDs

combined_list$Dolphin.ID1 <- affil_females$Dolphin.ID[match(combined_list$ID1, affil_females$Combined.ID)]
combined_list$Dolphin.ID2 <- affil_females$Dolphin.ID[match(combined_list$ID2, affil_females$Combined.ID)]

combined_pairs <- merge_pairs(combined_list, pairwise, 
                              xID1 = "Dolphin.ID1", xID2 = "Dolphin.ID2", 
                              yID1 = "Dolphin1", yID2 = "Dolphin2", 
                              all.x = TRUE, all.y = FALSE)

combined_pairs <- combined_pairs[which(combined_pairs$weight!= 0),]
combined_pairs <- reduce_pairs(combined_pairs, ID1 = "ID1", ID2 = "ID2")

#recreate network from data frame

network <- graph_from_data_frame(combined_pairs[,c("ID1", "ID2", "weight", "DyadML", "Closekin")], 
                                 directed = FALSE)

kin_graph <- subgraph_from_edges(network, eids = E(network)[E(network)$Closekin == "Y"])

#Read in individual combined ID strength calculation 
tot_str <- read.csv("Outputs/total_strength.csv")

#Calculating kin strength
str_kin <- strength(kin_graph)
str_kin_df <- data.frame(name = V(kin_graph)$name, strength = str_kin)

#Merge strength calculations
affil_females <- merge(affil_females, str_kin_df, by.x = "Combined.ID", by.y = "name", all.x = TRUE)
affil_females <- affil_females %>% rename(kin.strength = strength)

affil_females <- merge(affil_females, tot_str, by.x = "Combined.ID", by.y = "name", all.x = TRUE)
affil_females <- affil_females %>% rename(tot.strength = strength) 
affil_females <- affil_females %>% filter(!is.na(kin.strength))

#Normalize by kin strength divided by total strength
affil_females$norm.kin <- affil_females$kin.strength / affil_females$tot.strength

#Filter out Repro category unknown
affil_females <- affil_females %>% filter(ReproCat != "unknown")

#write.csv(affil_females, "Outputs/affil_kin_strength.csv", row.names = FALSE)

#Plot Raw Close Kin Strength 
kinmeans_raw <- affil_females %>%
  group_by(ReproCat) %>%
  summarise(mean_value = mean(kin.strength, na.rm = TRUE))

ggplot(affil_females, aes(x = ReproCat, y = kin.strength, fill = ReproCat)) +
  geom_boxplot() +
  geom_segment(data = kinmeans_raw, aes(x = as.numeric(ReproCat) - 0.2, xend = as.numeric(ReproCat) + 0.2,
                                     y = mean_value, yend = mean_value), color = "black", size = 1) +
  geom_text(data = kinmeans_raw, aes(x = ReproCat, y = mean_value, label = round(mean_value, 2)),
            vjust = -0.5, color = "black") +
  labs(x = "Repro State", y = "Raw Close Kin Strength") +
  theme_minimal() +
  theme(legend.position = "none")

#Plot Normalized close kin strength
kinmeans_norm <- affil_females %>%
  group_by(ReproCat) %>%
  summarise(mean_value = mean(norm.kin, na.rm = TRUE))

ggplot(affil_females, aes(x = ReproCat, y = norm.kin, fill = ReproCat)) +
  geom_boxplot() +
  geom_segment(data = kinmeans_norm, aes(x = as.numeric(ReproCat) - 0.2, xend = as.numeric(ReproCat) + 0.2,
                                        y = mean_value, yend = mean_value), color = "black", size = 1) +
  geom_text(data = kinmeans_norm, aes(x = ReproCat, y = mean_value, label = round(mean_value, 2)),
            vjust = -0.5, color = "black") +
  labs(x = "Repro State", y = "Normalized Close Kin Strength") +
  theme_minimal() +
  theme(legend.position = "none")

#Plot Raw Close Kin Strength with age class
agemeans_raw <- affil_females %>%
  group_by(AgeCat) %>%
  summarise(mean_value = mean(kin.strength, na.rm = TRUE))

ggplot(affil_females, aes(x = AgeCat, y = kin.strength, fill = AgeCat)) +
  geom_boxplot() +
  geom_segment(data = agemeans_raw, aes(x = as.numeric(AgeCat) - 0.2, xend = as.numeric(AgeCat) + 0.2,
                                        y = mean_value, yend = mean_value), color = "black", size = 1) +
  geom_text(data = agemeans_raw, aes(x = AgeCat, y = mean_value, label = round(mean_value, 2)),
            vjust = -0.5, color = "black") +
  labs(x = "Age Class", y = "Raw Close Kin Strength") +
  theme_minimal() +
  theme(legend.position = "none")

#Plot Normalized close kin strength with age class
agemeans_norm <- affil_females %>%
  group_by(AgeCat) %>%
  summarise(mean_value = mean(norm.kin, na.rm = TRUE))

ggplot(affil_females, aes(x = AgeCat, y = norm.kin, fill = AgeCat)) +
  geom_boxplot() +
  geom_segment(data = agemeans_norm, aes(x = as.numeric(AgeCat) - 0.2, xend = as.numeric(AgeCat) + 0.2,
                                         y = mean_value, yend = mean_value), color = "black", size = 1) +
  geom_text(data = agemeans_norm, aes(x = AgeCat, y = mean_value, label = round(mean_value, 2)),
            vjust = -0.5, color = "black") +
  labs(x = "Age Class", y = "Normalized Close Kin Strength") +
  theme_minimal() +
  theme(legend.position = "none")
