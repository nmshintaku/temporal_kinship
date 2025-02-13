#03 Statistical Analysis 

library(readr)
library(car)
library(SocGen)
library(igraph)
library(dplyr)
library(tidyr)
library(ggplot2)
library(bisonR)
library(brms)

#age <- read.csv("Outputs/affil_age_strength.csv")
#kin <- read.csv("Outputs/affil_kin_strength.csv")
#repro <- read.csv("Outputs/affil_repro_strength.csv")

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
sightings <- sightings %>%
  mutate(Repro.ID = sub("^([^.]+\\.[^.]+).*", "\\1", Combined.ID)) %>%
  mutate(Age.ID = paste(sub("\\..*", "", Combined.ID), sub(".*\\.(\\w+)$", "\\1", Combined.ID), sep = "."))
#Adding separate columns for Repro.ID and Age.ID
affil_females <- affil_females %>%
  mutate(Repro.ID = sub("^([^.]+\\.[^.]+).*", "\\1", Combined.ID)) %>%
  mutate(Age.ID = paste(sub("\\..*", "", Combined.ID), sub(".*\\.(\\w+)$", "\\1", Combined.ID), sep = "."))

affil_sightings <- sightings[which(sightings$Combined.ID %in% affil_females$Combined.ID), ]

#Filter out repro sightings >= 5 from affil_sightings
filtered_sightings <- affil_sightings %>% filter(ReproSightings >= 5)
affil_females <- affil_females %>%
  filter(Combined.ID %in% filtered_sightings$Combined.ID)

# mask the data so that association rates are only estimated in the timeframe where both members are alive
affil_mask <- schedulize(affil_females,
                         id = "Repro.ID",
                         start = "entrydate",
                         end = "departdate",
                         dates = dates,
                         format = "mask")

#Creating ReproCat column
extract_state <- function(Repro.ID) {
  return(sub(".*\\.", "", Repro.ID))
}

affil_females$ReproCat <- sapply(affil_females$Repro.ID, extract_state)

#Calculating SRI
masked_network <- simple_ratio(sightings = filtered_sightings,
                               group_variable = "Observation.ID",
                               dates = "Observation.Date",
                               IDs = "Repro.ID",
                               mask = affil_mask)

# For each female reprocat combination, sum up per reprocat and per relatedness cat 

# Create a network object

masked_network[is.nan(masked_network)] <- 0
masked_network[is.na(masked_network)] <- 0

combined_list <- mat2dat(masked_network, value.name = "weight")

#Set kinship as the edge attribute 

#create list of all pairwise combinations of combined IDs

combined_list$Dolphin.ID1 <- affil_females$Dolphin.ID[match(combined_list$ID1, affil_females$Repro.ID)]
combined_list$Dolphin.ID2 <- affil_females$Dolphin.ID[match(combined_list$ID2, affil_females$Repro.ID)]

combined_pairs <- merge_pairs(combined_list, pairwise, 
                              xID1 = "Dolphin.ID1", xID2 = "Dolphin.ID2", 
                              yID1 = "Dolphin1", yID2 = "Dolphin2", 
                              all.x = TRUE, all.y = FALSE)

combined_pairs <- combined_pairs[which(combined_pairs$weight!= 0),]
combined_pairs <- reduce_pairs(combined_pairs, ID1 = "ID1", ID2 = "ID2")

#recreate network from data frame

network <- graph_from_data_frame(combined_pairs[,c("ID1", "ID2", "weight", "DyadML", "Closekin")], 
                                 directed = FALSE)

synchrony_graph <- subgraph_from_edges(network, eids = E(network)[E(network)$Closekin == "Y"])

#################
#BISoN Regression
#################

#Run SRI with association = x to get count of times seen together
pair_count <- simple_ratio(sightings = filtered_sightings,
                               group_variable = "Observation.ID",
                               dates = "Observation.Date",
                               IDs = "Repro.ID",
                               mask = affil_mask,
                               assocInd = "X")

pair_count[is.nan(pair_count)] <- 0
pair_count[is.na(pair_count)] <- 0

pair_count_list <- mat2dat(pair_count, value.name = "pair.count")

combined_pairs <- merge_pairs(combined_pairs, pair_count_list, 
                              xID1 = "ID1", xID2 = "ID2", 
                              yID1 = "ID1", yID2 = "ID2", 
                              all.x = TRUE, all.y = FALSE)

#Run SRI with association = xy to get total count of joined sightings
total_count <- simple_ratio(sightings = filtered_sightings,
                           group_variable = "Observation.ID",
                           dates = "Observation.Date",
                           IDs = "Repro.ID",
                           mask = affil_mask,
                           assocInd = "XY")

total_count[is.nan(total_count)] <- 0
total_count[is.na(total_count)] <- 0

total_count_list <- mat2dat(total_count, value.name = "total.count")

combined_pairs <- merge_pairs(combined_pairs, total_count_list, 
                              xID1 = "ID1", xID2 = "ID2", 
                              yID1 = "ID1", yID2 = "ID2", 
                              all.x = TRUE, all.y = FALSE)

#Create edge model 
fit_edge <- bison_model(
  (pair.count | total.count) ~ dyad(ID1, ID2),
  data = combined_pairs,
  model_type = "count")

summary(fit_edge)

#Run dyadic regression
dyadic <- bison_brm(
  bison(edge_weight(ID1, ID2)) ~ DyadML + (1|mm(ID1, ID2)),
  fit_edge,
  combined_pairs,
  chains = 2,
  cores = 4
)


###########
##ANOVA####
###########

age <- read.csv("Outputs/affil_age_strength.csv")
kin <- read.csv("Outputs/affil_kin_strength.csv")
repro <- read.csv("Outputs/affil_repro_strength.csv")


#Repro ANOVA
repro_anova <- aov(norm.repro.indiv ~ ReproCat, data = repro)
summary(repro_anova)
plot(repro_anova)
repro_tukey <- TukeyHSD(repro_anova)
print(repro_tukey)
plot(repro_tukey)

#Close Kin ANOVA
kin_anova <- aov(DyadML ~ Closekin, data = combined_pairs)
summary(kin_anova)
plot(kin_anova)
kin_tukey <- TukeyHSD(kin_anova)
print(kin)

#Need to use combined_pairs (with combined ID) with an age category 
kin_anova_interaction <- aov(DyadML ~ Closekin * AgeCat, data = )

#Age ANOVA
age_anova <- aov(norm.age ~ AgeCat, data = age)
summary(age_anova)
plot(age_anova)
age_tukey <- TukeyHSD(age_anova)
print(age_tukey)
