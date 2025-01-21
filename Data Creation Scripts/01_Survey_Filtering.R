#01 Survey Dolphin Filtering 

library(tidyverse)
library(SocGen) ## remotes::install_github("vjf2/SocGen")
library(readxl)
vignette("dplyr")

#####################
#Survey Dolphin File
#####################

#Using survey dolphin file, filtering out low for dolphin certainty and true for beyond 10m
survey <- read.csv("Raw Data/SurveyDolphin_20250116.csv")

#Should this be & or | symbol?
#should be 204,233 observations?
filtered_survey <- survey %>% 
  filter((Dolphin.ID.Certainty == "HIGH" & Dolphin.ID.Certainty == "") |
           Beyond.10.Meters.Ind == "False")

#filtered_survey <- survey %>% filter(Dolphin.ID.Certainty == "HIGH" | Dolphin.ID.Certainty == "")
#filtered_survey <- survey %>% filter(Beyond.10.Meters.Ind == "False")

#Filter 1 sighting per day per individual using last sighting
filtered_survey <- filtered_survey %>%
  group_by(Dolphin.ID, Observation.Date) %>%
  arrange(desc(Observation.End.Time)) %>%
  slice(1) %>%
  ungroup()

#Keep the columns we need
filtered_survey <- filtered_survey %>%
  select(Observation.ID, Observation.Date, Dolphin.ID)

#write.csv(filtered_survey, "Outputs/filtered_survey.csv", row.names = FALSE)




