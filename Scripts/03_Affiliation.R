####################################
#### Classify and Model Affiliations
####################################

library(SocGen)
library(sna)
library(car)

load("IntermediateData/affiliation_networks.RData")

# Convert association matrices to data frame and combine real and random data

masked_network[is.nan(masked_network)] <- 0

mn <- mat2dat(masked_network, "realHWI")

mr <- lapply(ai_affil_rand, mat2dat)

for (i in 1:length(mr)) {names(mr[[i]])[3] <- paste0("randHWI", i)}

mr <- c(list(mn), mr)

# Merge removes anything missing, individuals that didn't overlap have NaN and get removed in merge
fdata <- Reduce(function(x, y) merge(x, y, all.x = T), mr) # about 2 minutes

end_row <- ncol(fdata)

# Pull out anything above 0.975 as affiliation
fdata$quantile975 <- apply(fdata, 1, function(x) quantile(as.numeric(x[4:end_row]), 0.975, na.rm = TRUE))

fdata$affiliation975 <- ifelse(fdata$realHWI > fdata$quantile975, TRUE, FALSE)

# Calculate median and mean random HWI as expected values
fdata$medianHWI <- apply(fdata[, 4:end_row], 1, median, na.rm = TRUE)
fdata$meanHWI <- apply(fdata[, 4:end_row], 1, mean, na.rm = TRUE)

# Remove random values
fdata <- fdata[, grep("rand", names(fdata), invert = TRUE)]

# Add in pairwise covariates
pair_covar <- read.csv("SharedData/pairwise_covariates.csv")

fdata <- merge_pairs(fdata, pair_covar, xID1 = "ID1", xID2 = "ID2")

fdata <- fdata[!duplicated(fdata), ]

# write.csv(fdata, file = "IntermediateData/affiliations.csv")

# Generate matrices for MRQAP
moddata <- reduce_pairs(fdata, "ID1", "ID2")

# Function to convert vectors to matrices
matricize <- function(x) {
  g <- igraph::graph.data.frame(x, directed = FALSE)
  m <- igraph::get.adjacency(g, attr = names(x)[3], sparse = FALSE)
  mat <- m[sort(rownames(m)), sort(colnames(m))]
  diag(mat) <- 1
  return(mat)
}

model_matrices <- moddata

model_matrices$ta <- ifelse(model_matrices$affiliation975 == TRUE, 1, 0)

# Set response to missing if kin status is unknown or if yrs_overlap is less than 1
model_matrices$ta <- ifelse(model_matrices$kintype == "close unknown", NA, model_matrices$ta)
model_matrices$ta <- ifelse(model_matrices$yrsoverlap <= 1, NA, model_matrices$ta)

# Format variables
model_matrices$m2 <- ifelse(model_matrices$kintype == "maternal", 1, 0)
model_matrices$p3 <- ifelse(model_matrices$kintype == "paternal", 1, 0)
model_matrices$d4 <- ifelse(model_matrices$kintype == "distant", 1, 0)

ta <- matricize(model_matrices[, c("ID1", "ID2", "ta")])
m2 <- matricize(model_matrices[, c("ID1", "ID2", "m2")])
p3 <- matricize(model_matrices[, c("ID1", "ID2", "p3")])
d4 <- matricize(model_matrices[, c("ID1", "ID2", "d4")])

hr <- matricize(model_matrices[, c("ID1", "ID2", "HRO")])
yr <- matricize(model_matrices[, c("ID1", "ID2", "yrsoverlap")])
js <- matricize(model_matrices[, c("ID1", "ID2", "joint_sightings")])

# Run MRQAP
ta <- as.sociomatrix.sna(ta)

set.seed(286567440)

mod <- sna::netlogit(y = ta, x = list(m2, p3, d4, hr, yr, js),
                     mode = "graph",
                     test.statistic = "z-value",
                     intercept = TRUE, reps = 1000)

summary(mod)

# save(mod, file="IntermediateData/mrqap.RData") #run time ~5 min

# Create Table 1
results <- data.frame(Estimate = mod$coefficients,
                      OddsRatio = exp(mod$coefficients),
                      Std.Error = mod$se,
                      Zvalue = mod$tstat,
                      Zlower = apply(mod$dist, 2, function(x) quantile(x, 0.025)),
                      Zupper = apply(mod$dist, 2, function(x) quantile(x, 0.975)),
                      Pvalue = mod$pgreqabs)

results <- apply(results, 2, function(x) round(x, 3))

# Format Table
results <- as.data.frame(results)
results$Zvalue <- paste0(results$Zvalue, " (", results$Zlower, " - ", results$Zupper, ")")
results$Pvalue <- ifelse(results$Pvalue < 0.001, "<0.001", results$Pvalue)
results[1, "Pvalue"] <- "-"

rownames(results) <- c("Intercept",
                       "Close Matrilineal Kin",
                       "Close Non-Matrilineal Kin",
                       "Distant Kin",
                       "Home Range Overlap (BA)",
                       "Years Observed",
                       "Joint # Observations")

# write.csv(results[,c(1:4,7)], "TablesFigures/Table1.csv")

# Extract model matrix from mrqap model (just to confirm same model specification)

x = list(m2, p3, d4, hr, yr, js)
nx <- stackcount(x) + 1
y <- ta
n <- dim(y)[2]
g <- list(y)

g[[2]] <- matrix(1, n, n)

for (i in 1:(nx - 1)) g[[i + 1 + 1]] <- x[[i]]

gv <- lapply(g[1:length(g)], function(x) gvectorize(x, mode = "graph", diag = FALSE))

mm <- do.call("cbind", gv)

mm <- mm[complete.cases(mm), ]

mm <- as.data.frame(mm)

names(mm) <- c("response", "intercept", "maternal",
               "paternal", "distant", "HRO", "yrs", "js")

glmod <- glm(response ~ maternal + paternal + distant + HRO + yrs + js,
             family = "binomial",
             data = mm)

vcov(glmod)
car::vif(glmod)

summary(glmod)

# Save glmod output for Figure 2
# save(glmod, file="IntermediateData/standardformmod.RData")