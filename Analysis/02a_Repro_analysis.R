#02a Association Scratch script (Repro)

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

#Merge repro strength to affil females dataframe
affil_females <- merge(affil_females, str_repro_df, by.x = "Repro.ID", by.y = "name", all.x = TRUE)
affil_females <- affil_females %>% rename(repro.strength = strength)

#Setting vertex attributes to each reproductive category
repro_network <- set_vertex_attr(repro_network, "ReproCat",
                                 value = sightings$ReproCat[match(V(repro_network)$name, sightings$Repro.ID)])


###############
#Repro Analysis
###############

#Average raw value per repro category across all females
lact <- affil_females %>%
  filter(grepl("\\.lact$", Repro.ID)) %>%
  summarise(lact = mean(repro.strength))
  
cyc <- affil_females %>%
  filter(grepl("\\.cyc$", Repro.ID)) %>%
  summarise(cyc = mean(repro.strength))

preg <- affil_females %>%
  filter(grepl("\\.preg$", Repro.ID)) %>%
  summarise(preg = mean(repro.strength))

juvenile <- affil_females %>%
  filter(grepl("\\.juvenile$", Repro.ID)) %>%
  summarise(juvenile = mean(repro.strength))

avg_repro_strength <- bind_cols(cyc, juvenile, lact, preg)

#STEP 1 Normalize avg value of each repro category 
#Dividing repro.strength of each individual by the average of the respective repro category 

extract_state <- function(Repro.ID) {
  return(sub(".*\\.", "", Repro.ID))
}

affil_females$ReproCat <- sapply(affil_females$Repro.ID, extract_state)

affil_females$norm_repro <- mapply(function(repro.strength, ReproCat) {
  avg_value <- avg_repro_strength[[ReproCat]]
  return(repro.strength / avg_value)
}, affil_females$repro.strength, affil_females$ReproCat)

#Code not working to get rid of character 0
affil_females$norm_repro[affil_females$norm_repro == "character(0)"] <- NA

#filter out unknown ReproCat
no_unk_affil <- subset(affil_females, ReproCat != "unknown")
no_unk_affil$norm_repro <- as.numeric(no_unk_affil$norm_repro) 


#STEP 2 Normalizing by individual level

#Calculating total average strength per individual across all states
no_unk_affil <- no_unk_affil %>%
  group_by(Dolphin.ID) %>%
  mutate(avg.tot.str = mean(repro.strength),
         norm.repro.indiv = repro.strength / avg.tot.str)
#No Nan or inf
clean_affil <- no_unk_affil %>%
  filter(!is.nan(norm.repro.indiv) & !is.infinite(norm.repro.indiv))

#Plot repro strength normalized to individual with average shown across each state
means_indiv <- clean_affil %>%
  group_by(ReproCat) %>%
  summarise(mean_value = mean(norm.repro.indiv, na.rm = TRUE))

ggplot(clean_affil, aes(x = ReproCat, y = norm.repro.indiv, fill = ReproCat)) +
  geom_boxplot() +
  geom_segment(data = means_indiv, aes(x = as.numeric(ReproCat) - 0.2, xend = as.numeric(ReproCat) + 0.2,
                                 y = mean_value, yend = mean_value), color = "black", size = 1) +
  geom_text(data = means_indiv, aes(x = ReproCat, y = mean_value, label = round(mean_value, 2)),
            vjust = -0.5, color = "black") +
  labs(x = "Repro State", y = "Normalized Strength to Individual") +
  theme_minimal() +
  theme(legend.position = "none")


#Repro strength normalized strength to the category level
means_norm <- clean_affil %>%
  group_by(ReproCat) %>%
  summarise(mean_value = mean(norm_repro, na.rm = TRUE))

ggplot(clean_affil, aes(x = ReproCat, y = norm_repro, fill = ReproCat)) +
  geom_boxplot() +
  geom_segment(data = means_norm, aes(x = as.numeric(ReproCat) - 0.2, xend = as.numeric(ReproCat) + 0.2,
                                       y = mean_value, yend = mean_value), color = "black", size = 1) +
  geom_text(data = means_norm, aes(x = ReproCat, y = mean_value, label = round(mean_value, 2)),
            vjust = -0.5, color = "black") +
  labs(x = "Repro State", y = "Normalized Strength to Category") +
  theme_minimal() +
  theme(legend.position = "none")

#Raw Repro Strength
means_raw <- clean_affil %>%
  group_by(ReproCat) %>%
  summarise(mean_value = mean(repro.strength, na.rm = TRUE))

ggplot(clean_affil, aes(x = ReproCat, y = repro.strength, fill = ReproCat)) +
  geom_boxplot() +
  geom_segment(data = means_raw, aes(x = as.numeric(ReproCat) - 0.2, xend = as.numeric(ReproCat) + 0.2,
                                      y = mean_value, yend = mean_value), color = "black", size = 1) +
  geom_text(data = means_raw, aes(x = ReproCat, y = mean_value, label = round(mean_value, 2)),
            vjust = -0.5, color = "black") +
  labs(x = "Repro State", y = "Raw Repro Strength") +
  theme_minimal() +
  theme(legend.position = "none")


#Scratch

filtered_sightings <- affil_sightings %>% filter(ReproSightings >= 5)
affil_females <- affil_females %>%
  filter(Combined.ID %in% filtered_sightings$Combined.ID)

affil_females <- affil_females %>%
  group_by(Dolphin.ID) %>%
  mutate(avg.tot.repro.str = mean(repro.strength),
         norm.repro.indiv = repro.strength / avg.tot.repro.str)

affil_females[] <- lapply(affil_females, function(x) {
  if (is.list(x)) {
    sapply(x, toString)
  } else {
    x
  }
})

write.csv(affil_females, "Outputs/affil_repro_strength.csv", row.names = FALSE)

