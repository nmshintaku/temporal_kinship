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
  filter(age > 4)

#Calculate total number of sightings 
totalsightings <- adultF %>%
  group_by(Dolphin.ID) %>%
  summarise(Sightings = n_distinct(Observation.Date)) %>%
  ungroup()
#Now filtering for minimum of 10 sightings
filtered_sightings <- totalsightings %>%
  filter(Sightings >= 10)
#write.csv(filtered_sightings, "Outputs/filtered_sightings.csv", row.names = FALSE)

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

individual_covariates <- individual_covariates %>%
  mutate(deathdate = if_else(is.na(deathdate), as.Date("2024-11-30"), deathdate))

#OLD STEPS for reference only 
#Filtering out individual's last observation date of the year 2024
#obs2024 <- filtered_survey %>%
 # filter(year(Observation.Date) == 2024) %>%
  #group_by(Dolphin.ID) %>%
  #summarise(LastObs = max(Observation.Date)) %>%
  #ungroup()
#if depart/death date is na, put 2024 last observation date
#individual_covariates <- individual_covariates %>% 
 # left_join(obs2024, by = "Dolphin.ID") %>%
  #mutate(deathdate = if_else(is.na(deathdate), LastObs, deathdate)) %>%
  #select(-LastObs)

individual_covariates <- individual_covariates %>% rename(departdate = deathdate)
individual_covariates <- individual_covariates %>% rename(entrydate = Turned4Date)

write.csv(individual_covariates, "Outputs/individual_covariates.csv", row.names = FALSE)


##################
#Sightings
##################

adultF$Observation.Date <- as.Date(adultF$Observation.Date)

#Creating age class
adultF <- adultF %>%
  mutate(AgeClass = case_when(
    age >= 4 & age < 10 ~ "juvenile",
    age >= 10 & age <= 30 ~ "adult",
    age > 30 ~ "elder",
    TRUE ~ NA_character_
  ))

#Creating reproductive categories
adultF <- adultF %>%
  mutate(ReproCat = case_when(
    pregnant == 1 ~ "preg",
    lactating == 1 ~ "lact",
    cycling == 1 ~ "cyc",
    TRUE ~ NA_character_
  )) %>%
  mutate(ReproCat = case_when(
    !is.na(ReproCat) ~ ReproCat,
    is.na(ReproCat) & AgeClass == "juvenile" ~ "juvenile",
    TRUE ~ "unknown"
  ))
  
sightings <- adultF %>%
  select(Dolphin.ID, Observation.Date) %>%
  left_join(filtered_survey %>% select(Dolphin.ID, Observation.Date, Observation.ID), by = c("Dolphin.ID", "Observation.Date")) %>%
  left_join(adultF %>% select(Dolphin.ID, Observation.Date, AgeClass, ReproCat),
            by = c("Dolphin.ID", "Observation.Date")) %>%
  mutate(Combined.ID = paste(Dolphin.ID, ReproCat, AgeClass, sep = ".")) %>%
  group_by(ReproCat) %>%
  mutate(ReproSightings = n()) %>%
  ungroup()

#Making sure to filter individuals with minimum 10 total sightings
sightings <- sightings %>%
  filter(Dolphin.ID %in% filtered_sightings$Dolphin.ID)

#calculate repro sightings per individual
rs_tab <- table(sightings$Combined.ID)
sightings$ReproSightings <- rs_tab[match(sightings$Combined.ID, names(rs_tab))]

#Keep repro sightings greater than or equal to 5
sightings <- sightings %>% filter(ReproSightings >= 5)

#Adding separate columns for Repro.ID and Age.ID
sightings <- sightings %>%
  mutate(Repro.ID = sub("^([^.]+\\.[^.]+).*", "\\1", Combined.ID)) %>%
  mutate(Age.ID = paste(sub("\\..*", "", Combined.ID), sub(".*\\.(\\w+)$", "\\1", Combined.ID), sep = "."))

write.csv(sightings, "Outputs/sightings.csv", row.names = FALSE)
#write.csv(adultF, "outputs/adultF_filtered.csv", row.names = FALSE)




