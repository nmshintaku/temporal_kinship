#02 Repro Assignment

library(tidyverse)
library(SocGen) ## remotes::install_github("vjf2/SocGen")
library(readxl)

#####################
#Life History File
#####################

lh <- read.csv("Raw Data/LifeHistory_20250116.csv")
filtered_survey <- read.csv("outputs/filtered_survey.csv")

lh$Birth.Date <- as.Date(lh$Birth.Date)
lh$Death.Date <- as.Date(lh$Death.Date)
filtered_survey$Observation.Date <- as.Date(filtered_survey$Observation.Date)

#Filter out males and unknown sex individuals 
allF <- lh$Dolphin.ID[which(lh$Sex == "FEMALE")]

allF <- filtered_survey[which(filtered_survey$Dolphin.ID %in% allF), c("Dolphin.ID", "Observation.Date")]

allF <- allF[!duplicated(allF), ]

# Add birthdate into allF
allF$birthdate <- lh$Birth.Date[match(
  allF$Dolphin.ID,
  lh$Dolphin.ID
)]

####################
# Define reproductive state
####################

calves <- lh[which(lh$Mother.ID != ""), ]

calves$pregnantstart <- calves$Birth.Date - 365
calves$cycstart6m <- calves$pregnantstart - (365 / 2)

calves$cycresume <- calves$Birth.Date + (365 * 1.8)

calves$cycresumeD <- ifelse(!is.na(calves$Death.Date) & calves$Death.Date < calves$cycresume, calves$Death.Date, calves$cycresume)

calves$cycresumeD <- as.Date(calves$cycresumeD, origin = "1970-01-01")

calves <- calves[, c("Dolphin.ID", "Dolphin.Name", "Mother.ID", "Birth.Date", "Death.Date", "pregnantstart", "cycstart6m", "cycresume", "cycresumeD")]

dates <- sort(unique(allF$Observation.Date))

pregnant <- schedulize(
  data = calves, id = "Dolphin.ID",
  start = "pregnantstart", end = "Birth.Date",
  dates = dates,
  format = c("mask")
)

colnames(pregnant) <- as.character(sort(unique(allF$Observation.Date)))
preg <- mat2dat(pregnant, "pregnant")
names(preg)[1] <- "pregID"

preg$mom <- calves$Mother.ID[match(preg$pregID, calves$Dolphin.ID)]


lactating <- schedulize(
  data = calves, id = "Dolphin.ID",
  start = "Birth.Date", end = "cycresumeD",
  dates = dates,
  format = c("mask")
)

colnames(lactating) <- as.character(sort(unique(allF$Observation.Date)))
lact <- mat2dat(lactating, "lactating")
names(lact)[1] <- "lactID"

lact$mom <- calves$Mother.ID[match(lact$lactID, calves$Dolphin.ID)]


defcyc <- schedulize(
  data = calves, id = "Dolphin.ID",
  start = "cycstart6m", end = "pregnantstart",
  dates = dates,
  format = c("mask")
)

colnames(defcyc) <- as.character(sort(unique(allF$Observation.Date)))
cyc <- mat2dat(defcyc, "cycling")
names(cyc)[1] <- "cycID"

cyc$mom <- calves$Mother.ID[match(cyc$cycID, calves$Dolphin.ID)]

# Merge into allF, then count adult males and merge total group size
allF <- merge(allF, preg,
              by.x = c("Dolphin.ID", "Observation.Date"),
              by.y = c("mom", "ID2"), all.x = TRUE
)

allF <- merge(allF, lact,
              by.x = c("Dolphin.ID", "Observation.Date"),
              by.y = c("mom", "ID2"), all.x = TRUE
)

allF <- merge(allF, cyc,
              by.x = c("Dolphin.ID", "Observation.Date"),
              by.y = c("mom", "ID2"), all.x = TRUE
)

# Add age

allF$age <- as.numeric(allF$Observation.Date - allF$birthdate) / 365.25

allF$mature <- ifelse(allF$age >= 11, "yes", "no")
allF$mature <- ifelse(!is.na(allF$cycling), "yes", allF$mature)
allF$mature <- ifelse(!is.na(allF$pregnant), "yes", allF$mature)
allF$mature <- ifelse(!is.na(allF$lactating), "yes", allF$mature)

allF$possiblycyc <- ifelse(is.na(allF$lactating) &
                             is.na(allF$pregnant) &
                             allF$mature == "yes", 1, NA)

# If known pregnancy overlaps with cycling of next calf, classify as only pregnant 

allF[is.na(allF)] <- 0

allF$cycling <- ifelse((allF$cycling==1 & allF$pregnant == 1), 
                       0, 
                       allF$cycling)

# If seen the same day as birth, classify as lactating

allF$pregnant <- ifelse((allF$pregnant==1 & allF$lactating == 1), 
                        0, 
                        allF$pregnant)


## End reproductive status determination here

write.csv(allF, "Outputs/reproductive_status_females_20250121.csv", row.names = FALSE)
