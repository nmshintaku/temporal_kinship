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
library(multcompView)

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

#Creating ReproCat column
extract_state <- function(Repro.ID) {
  return(sub(".*\\.", "", Repro.ID))
}

affil_females$ReproCat <- sapply(affil_females$Repro.ID, extract_state)

affil_females <- affil_females %>% filter(ReproCat != "unknown")

affil_sightings <- sightings[which(sightings$Combined.ID %in% affil_females$Combined.ID), ]
affil_sightings <- affil_sightings %>% filter(ReproCat != "unknown")

# mask the data so that association rates are only estimated in the timeframe where both members are alive
affil_mask <- schedulize(affil_females,
                         id = "Repro.ID",
                         start = "entrydate",
                         end = "departdate",
                         dates = dates,
                         format = "mask")

#Calculating SRI
masked_network <- simple_ratio(sightings = affil_sightings,
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

network <- set_vertex_attr(network, "ReproCat",
                           value = sightings$ReproCat[match(V(network)$name, sightings$Repro.ID)])

synchrony_graph <- subgraph_from_edges(network, eids = E(network)[E(network)$Closekin == "Y"])

#plot close kin + lact networks
lact_vertex <- V(network)[ReproCat == "lact"]
lact_graph <- induced.subgraph(network, vids = lact_vertex)
kin_edges <- E(lact_graph)[Closekin == "Y"]
kin_lact <- subgraph_from_edges(lact_graph, eids = kin_edges, delete.vertices = TRUE)
plot(kin_lact, edge.width = E(kin_lact)$weight)

#################
#BISoN Regression
#################

#Run SRI with association = x to get count of times seen together
pair_count <- simple_ratio(sightings = affil_sightings,
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
total_count <- simple_ratio(sightings = affil_sightings,
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

### Testing a subset of the data on Bison Model
#Creating subset of combined pairs to try on the Bison model 
subset_pairs <- combined_pairs[sample(nrow(combined_pairs), 20), ]
subset_edge <- bison_model(
  (pair.count | total.count) ~ dyad(ID1, ID2),
  data = subset_pairs,
  model_type = "count")

subset_dyadic <- bison_brm(
  bison(edge_weight(ID1, ID2)) ~ DyadML,
  fit_edge,
  subset_pairs,
  num_draws = 5,
  refresh = 0
)

### Bison model on full dataset
#Create edge model 
fit_edge <- bison_model(
  (pair.count | total.count) ~ dyad(ID1, ID2),
  data = combined_pairs,
  model_type = "count")

summary(fit_edge)

#Run dyadic regression
dyadic <- bison_brm(
  bison(edge_weight(ID1, ID2)) ~ DyadML,
  fit_edge,
  combined_pairs,
  num_draws = 5,
  refresh = 0
)

###Running bison simulation model 
sim_data <- simulate_bison_model("binary", aggregated = TRUE)
df <- sim_data$df_sim

priors <- get_default_priors("binary_conjugate")
prior_check(priors, "binary_conjugate")
priors$edge = "normal(-1, 2.5)"
fit_edge_sim <- bison_model(
  (event | duration) ~ dyad(node_1_id, node_2_id),
  data = df,
  model_type = "count"
)

df_dyadic <- df %>%
  distinct(node_1_id, node_2_id, age_diff)

fit_edge_sim <- bison_brm(
  bison(edge_weight(node_1_id, node_2_id)) ~ age_diff,
  fit_edge_sim,
  df, 
  num_draws = 5,
  refresh = 0
)


###########
##ANOVA####
###########

age <- read.csv("Outputs/affil_age_strength.csv")
kin <- read.csv("Outputs/affil_kin_strength.csv")
repro <- read.csv("Outputs/affil_repro_strength.csv")

#Filter out Repro category unknown
repro <- repro %>% filter(ReproCat != "unknown")
repro <- repro %>%
  filter(!is.na(norm.repro.indiv) & !is.infinite(norm.repro.indiv))

#Repro ANOVA
repro_anova <- aov(norm.repro.indiv ~ ReproCat, data = repro)
summary(repro_anova)
plot(repro_anova)
repro_tukey <- TukeyHSD(repro_anova)

print(repro_tukey)
plot(repro_tukey)

means <- repro %>%
  group_by(ReproCat) %>%
  summarize(mean_value = mean(norm.repro.indiv))
#To get group combinations means, add the individual group means and divide by number of groups 
#Juvenile - cyc = (1.5465586+0.9826508)/2 = 1.264605
tukey_df <- as.data.frame(repro_tukey$ReproCat)
tukey_df$comparison <- rownames(tukey_df)
tukey_df <- tukey_df %>%
  rename(
    p_value = `p adj`
  )

ggplot(tukey_df, aes(x = comparison, y = diff)) +
  geom_boxplot() +
  geom_text(data = tukey_df, aes(x = comparison, y = max(repro$norm.repro.indiv) + 0.5, label = round(p_value, 3)), color = "red") +
  labs(title = "Boxplots with Tukey HSD Results", x = "Reproductive Categories", y = "Values") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#Close Kin ANOVA
kin_anova <- aov(DyadML ~ Closekin, data = combined_pairs)
summary(kin_anova)
plot(kin_anova)
kin_tukey <- TukeyHSD(kin_anova)
print(kin)

#Creating ReproCat column
extract_age <- function(Age.ID) {
  return(sub(".*\\.", "", Age.ID))
}

kin$AgeCat <- sapply(kin$Age.ID, extract_age)

kin_age <- aov(norm.kin ~ AgeCat, data = kin)
summary(kin_age)

#Age ANOVA
age_anova <- aov(norm.age ~ AgeCat, data = age)
summary(age_anova)
plot(age_anova)
age_tukey <- TukeyHSD(age_anova)
print(age_tukey)

#Repro ANOVA, same state vs different state results
shared_state <- function(dolphin_id) {
  return(strsplit(dolphin_id, "\\.")[[1]][2])
}

# Apply the function to create new columns for states
repro_combined <- combined_pairs %>%
  mutate(State1 = sapply(ID1, shared_state),
         State2 = sapply(ID2, shared_state))

# Check if the states are the same and create the shared_state column
repro_combined <- repro_combined %>%
  mutate(shared_state = ifelse(State1 == State2, "Y", "N"))

# Drop the temporary state columns
repro_combined <- repro_combined %>%
  select(-State1, -State2)

repro_state <- aov(weight ~ shared_state, data = repro_combined)
summary(repro_state)

repro_state_interaction <- aov(weight ~ shared_state * Closekin, data = repro_combined)
summary(repro_state_interaction)
repro_results <- Anova(repro_state_interaction, type = "II")
p_value <- repro_results$`Pr(>F)`[1]

ggplot(repro_combined, aes(x = shared_state, y = weight)) +
  geom_boxplot(outlier.shape = 16, outlier.size = 2) +
  stat_summary(fun = mean, geom = "point", shape = 16, size = 3, fill = "black") +
  stat_summary(fun = mean, geom = "text", aes(label = round(..y.., 2)), vjust = -0.5) +
  labs(x = "Shared State", y = "Weight") +
  annotate("text", x = 1.5, y = max(repro_combined$weight), label = paste("p-value:", round(p_value, 4)), size = 5, hjust = 0.5) +
  theme_minimal()

#Reporting Statistics
#Number of unique individuals in affil females
length(unique(affil_females$Dolphin.ID)) #181
length(unique(affil_sightings$Dolphin.ID)) #181

#Average sighting per dolphin
average_sightings <- affil_females %>%
  distinct(Dolphin.ID, .keep_all = TRUE) %>%
  select(Dolphin.ID, Sightings)
mean(average_sightings$Sightings) #112
median(average_sightings$Sightings) #55
print(IQR(average_sightings$Sightings))
print(quantile(average_sightings$Sightings))
