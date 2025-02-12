#02a Association Scratch script (Age)

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

#Adding separate columns for Repro.ID and Age.ID
sightings <- sightings %>%
  mutate(Repro.ID = sub("^([^.]+\\.[^.]+).*", "\\1", Combined.ID)) %>%
  mutate(Age.ID = paste(sub("\\..*", "", Combined.ID), sub(".*\\.(\\w+)$", "\\1", Combined.ID), sep = "."))

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

#Adding separate columns for Repro.ID and Age.ID
affil_sightings <- affil_sightings %>%
  mutate(Repro.ID = sub("^([^.]+\\.[^.]+).*", "\\1", Combined.ID)) %>%
  mutate(Age.ID = paste(sub("\\..*", "", Combined.ID), sub(".*\\.(\\w+)$", "\\1", Combined.ID), sep = "."))

#Creating AgeCat column
extract_state <- function(Age.ID) {
  return(sub(".*\\.", "", Age.ID))
}

affil_females$AgeCat <- sapply(affil_females$Age.ID, extract_state)

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

#Read in individual combined ID strength calculation 
tot_str <- read.csv("Outputs/total_strength.csv")

affil_females <- merge(affil_females, tot_str, by.x = "Combined.ID", by.y = "name", all.x = TRUE)
affil_females <- affil_females %>% rename(tot.strength = strength) 

#Merge strength calculations
affil_females <- merge(affil_females, str_age_df, by.x = "Age.ID", by.y = "name", all.x = TRUE)
affil_females <- affil_females %>% rename(age.strength = strength)

#Setting vertex attributes to each age class
age_network <- set_vertex_attr(age_network, "AgeClass",
                               value = sightings$AgeClass[match(V(age_network)$name, sightings$Age.ID)])


###############
#Age Analysis
###############

#Filter out repro sightings >= 5 from affil_sightings
filtered_sightings <- affil_sightings %>% filter(ReproSightings >= 5)
affil_females <- affil_females %>%
  filter(Combined.ID %in% filtered_sightings$Combined.ID)

#Normalize by kin strength divided by total strength
affil_females$norm.age <- affil_females$age.strength / affil_females$tot.strength

#Calculating total average strength per individual across all age classes
affil_females <- affil_females %>%
  group_by(Dolphin.ID) %>%
  mutate(avg.tot.age.str = mean(age.strength),
         norm.age.indiv = age.strength / avg.tot.age.str)

#write.csv(affil_females, "Outputs/affil_age_strength.csv", row.names = FALSE)

#Plot Raw Age Strength 
agemeans_raw <- affil_females %>%
  group_by(AgeCat) %>%
  summarise(mean_value = mean(age.strength, na.rm = TRUE))

ggplot(affil_females, aes(x = AgeCat, y = age.strength, fill = AgeCat)) +
  geom_boxplot() +
  geom_segment(data = agemeans_raw, aes(x = as.numeric(AgeCat) - 0.2, xend = as.numeric(AgeCat) + 0.2,
                                        y = mean_value, yend = mean_value), color = "black", size = 1) +
  geom_text(data = agemeans_raw, aes(x = AgeCat, y = mean_value, label = round(mean_value, 2)),
            vjust = -0.5, color = "black") +
  labs(x = "Age State", y = "Raw Age Strength") +
  theme_minimal() +
  theme(legend.position = "none")

#Plot Normalized age strength (raw age strength / tot strength)
agemeans_norm <- affil_females %>%
  group_by(AgeCat) %>%
  summarise(mean_value = mean(norm.age, na.rm = TRUE))

ggplot(affil_females, aes(x = AgeCat, y = norm.age, fill = AgeCat)) +
  geom_boxplot() +
  geom_segment(data = agemeans_norm, aes(x = as.numeric(AgeCat) - 0.2, xend = as.numeric(AgeCat) + 0.2,
                                         y = mean_value, yend = mean_value), color = "black", size = 1) +
  geom_text(data = agemeans_norm, aes(x = AgeCat, y = mean_value, label = round(mean_value, 2)),
            vjust = -0.5, color = "black") +
  labs(x = "Age State", y = "Normalized Age Strength") +
  theme_minimal() +
  theme(legend.position = "none")

#Plot age strength normalized to individual with average shown across each state
agemeans_indiv <- affil_females %>%
  group_by(AgeCat) %>%
  summarise(mean_value = mean(norm.age.indiv, na.rm = TRUE))

ggplot(affil_females, aes(x = AgeCat, y = norm.age.indiv, fill = AgeCat)) +
  geom_boxplot() +
  geom_segment(data = agemeans_indiv, aes(x = as.numeric(AgeCat) - 0.2, xend = as.numeric(AgeCat) + 0.2,
                                       y = mean_value, yend = mean_value), color = "black", size = 1) +
  geom_text(data = agemeans_indiv, aes(x = AgeCat, y = mean_value, label = round(mean_value, 2)),
            vjust = -0.5, color = "black") +
  labs(x = "Age State", y = "Normalized Age Strength to Individual") +
  theme_minimal() +
  theme(legend.position = "none")

