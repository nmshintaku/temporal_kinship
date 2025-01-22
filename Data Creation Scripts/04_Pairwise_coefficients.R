#04 Pairwise Coefficients

library(dplyr)
library(tidyr)

filtered_sightings <- read.csv("Outputs/filtered_sightings.csv")
relatedness <- read.csv("Raw Data/RelatednessEstimates_2024.csv")

#Calculate total number of sightings 
totalsightings <- adultF %>%
  group_by(Dolphin.ID) %>%
  summarise(Sightings = n_distinct(Observation.Date)) %>%
  ungroup()
#Now filtering for minimum of 10 sightings
filtered_sightings <- totalsightings %>%
  filter(Sightings >= 10)
#Checking dolphin numbers
#dolphins_count <- sightings %>%
 # summarise(value = n_distinct(Dolphin.ID))
#384 female dolphins in the project

unique_dolphins <- filtered_sightings %>%
  distinct(Dolphin.ID)
unique_dolphins$Dolphin.ID <- as.character(unique_dolphins$Dolphin.ID)

#Create all possible combination pairs

#THIS CODE NOT WORKING
#combinations <- expand.grid(Dolphin1 = unique_dolphins$Dolphin.ID,
                            #Dolphin2 = unique_dolphins$Dolphin.ID) %>%
 # filter(Dolphin1 < Dolphin2) %>%
  #mutate(combination = paste(Dolphin1, Dolphin2, sep = "."))


combinations <- as.data.frame(t(combn(unique_dolphins$Dolphin.ID, 2)))
colnames(combinations) <- c("Dolphin1", "Dolphin2")
combinations <- combinations %>%
  rowwise() %>%
  mutate(IDPair = paste(sort(c(Dolphin1, Dolphin2)), collapse = "."))

dyad <- relatedness %>%
  mutate(IDPair = as.character(IDPair)) %>%
  rowwise() %>%
  mutate(IDPair = paste(sort(strsplit(IDPair, "\\.")[[1]]), collapse = ".")) %>%
  select(IDPair, DyadML)
  

combinations <- combinations %>%
  left_join(dyad, by = "IDPair")

write.csv(combinations, "Outputs/pairwise_covariates.csv", row.names = FALSE)


