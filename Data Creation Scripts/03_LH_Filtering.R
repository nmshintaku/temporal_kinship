#03 Life History Filtering

library(tidyverse)
library(lubridate)

repro <- read.csv("Outputs/reproductive_status_females_20250121.csv")
lh <- read.csv("Raw Data/LifeHistory_20250116.csv")
filtered_survey <- read.csv("outputs/filtered_survey.csv")

repro$birthdate <- as.Date(repro$birthdate)
lh$Death.Date <- as.Date(lh$Death.Date)
filtered_survey$Observation.Date <- as.Date(filtered_survey$Observation.Date)

# Add deathdate into repro
repro$deathdate <- lh$Death.Date[match(
  repro$Dolphin.ID,
  lh$Dolphin.ID
)]

#Filter out calves; individuals older than 4 yrs
adultF <- repro %>%
  filter(mature == "yes")

#Calculate total number of sightings 
totalsightings <- adultF %>%
  group_by(Dolphin.ID) %>%
  summarise(Sightings = n_distinct(Observation.Date)) %>%
  ungroup()
#Now filtering for minimum of 10 sightings
filtered_sightings <- totalsightings %>%
  filter(Sightings >= 10)

##################
#Individual Covars
##################

individual_covariates <- filtered_sightings %>%
  select(Dolphin.ID, Sightings)

#new dataframe for birth and death date
dates <- adultF %>%
  select(Dolphin.ID, birthdate, deathdate) %>%
  distinct()
#calculate date turning 4 years
dates <- dates %>%
  mutate(Turned4Date = birthdate %>% ymd() %>% add_with_rollback(years(4)))

individual_covariates <- individual_covariates %>% left_join(dates %>% select(Dolphin.ID, Turned4Date, deathdate), by = "Dolphin.ID")

#Filtering out last observation date of the year 2024
obs2024 <- filtered_survey %>%
  filter(year(Observation.Date) == 2024) %>%
  group_by(Dolphin.ID) %>%
  summarise(LastObs = max(Observation.Date)) %>%
  ungroup()
#if depart/death date is na, put 2024 last observation date
individual_covariates <- individual_covariates %>% 
  left_join(obs2024, by = "Dolphin.ID") %>%
  mutate(deathdate = if_else(is.na(deathdate), LastObs, deathdate)) %>%
  select(-LastObs)

individual_covariates <- individual_covariates %>% rename(departdate = deathdate)

#THERE ARE STILL NAs IN DEPART DATE
write.csv(individual_covariates, "Outputs/individual_covariates.csv", row.names = FALSE)
