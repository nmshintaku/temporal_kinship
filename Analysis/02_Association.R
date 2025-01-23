#######################################
#### Calculate Association Coefficients
#######################################

library(SocGen)
library(igraph)

indiv_covars <- read.csv("Outputs/individual_covariates.csv")
indiv_covars$entrydate <- as.Date(indiv_covars$entrydate)
indiv_covars$departdate <- as.Date(indiv_covars$departdate)

sightings <- read.csv("Outputs/sightings.csv")

dates <- sort(unique(sightings$Observation.Date))

#count repro sightings

rs_tab <- table(sightings$Combined.ID)
sightings$ReproSightings <- rs_tab[match(sightings$Combined.ID, names(rs_tab))]

# select just genotyped individuals

relatedness <- read.csv("Raw Data/RelatednessEstimates_2024.csv")
genotypes <- unique(c(relatedness$ID1, relatedness$ID2))

indiv_covars$genotyped <- ifelse(indiv_covars$Dolphin.ID %in% genotypes, "Y", "N")

statuses <- sightings[,c("Dolphin.ID", "Combined.ID")] |> unique()

affil_females <- merge(statuses, indiv_covars, by = "Dolphin.ID", all.x = TRUE)

affil_females <- affil_females[which(affil_females$genotyped == "Y"), ]

# mask the data so that association rates are only estimated in the timeframe where both members are alive
affil_mask <- schedulize(affil_females,
                         id = "Combined.ID",
                         start = "entrydate",
                         end = "departdate",
                         dates = dates,
                         format = "mask")

affil_sightings <- sightings[which(sightings$Combined.ID %in% affil_females$Combined.ID), ]

masked_network <- half_weight(sightings = affil_sightings,
                              group_variable = "Observation.ID",
                              dates = "Observation.Date",
                              IDs = "Combined.ID",
                              mask = affil_mask)

