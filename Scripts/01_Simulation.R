
###############################
#### Simulate Null Model Groups
###############################

library(parallel)
library(foreach)
library(doParallel)
library(doRNG)
library(SocGen) # remotes::install_github("vjf2/SocGen")
library(sp)

num_sim <- 1000 # set number of simulations to run

load("Data/home_ranges.RData")
load("Data/surveyeffort.RData")

gridrad <- home_ranges@grid@cellsize[1] / 2

# make schedule

indiv_covars <- read.csv("Data/individual_covariates.csv")
indiv_covars$entry <- as.Date(indiv_covars$entry)
indiv_covars$depart <- as.Date(indiv_covars$depart)

dates <- names(daily_search_effort)
d <- length(dates)

schedule <- schedulize(indiv_covars, id = "Dolphin.ID", start = "entry", end = "depart",
                       dates = dates, format = "sim")

sightings <- read.csv("Data/sightings.csv")

numdol <- c(table(sightings$Date))

# set up cluster for parallelization
# timing depends on number of cores available
# 7 threads - 1000 sims in 17 minutes

cl <- makeCluster(detectCores() - 1)
clusterEvalQ(cl, {library(sp);library(SocGen)})
clusterExport(cl, c("d", "daily_search_effort", "home_ranges", "schedule", "num_sim", "numdol", "gridrad"))
registerDoParallel(cl)
registerDoRNG(seed = 286567440)

starttime <- Sys.time()

nest_days <- foreach(i = seq_len(d), .errorhandling = 'pass') %dopar% {
  
  bound <- daily_search_effort[i, ]
  nd <- numdol[i]
  dailygrid <- home_ranges[bound, , drop = TRUE]
  probweights <- colSums(dailygrid@data, na.rm = TRUE)
  probweights <- probweights[names(probweights) %in% colnames(schedule)[schedule[i, ] == TRUE]]
  dc <- coordinates(dailygrid)
  dgdf <- dailygrid@data
  holder <- replicate(num_sim, fast_random_points(probweights = probweights,
                                                  nd = nd,
                                                  dc = dc,
                                                  dgdf = dgdf,
                                                  gridrad = gridrad),
                      simplify = FALSE)
  return(holder)
}

endtime <- Sys.time()

stopCluster(cl)

endtime - starttime # run time ~20 min on Intel(R) Core(TM) i7-8565U CPU @ 1.80GHz

# rearrange list so that days are listed within iteration
sim_surveys <- sapply(1:num_sim, function(i) lapply(nest_days, "[[", i), simplify = FALSE)

# group dolphins together using hierarchical distance-based clustering to get same average group size as in real data
groupperday <- table(sightings$Date[!duplicated(sightings$Observation.ID)])

# group_assign run time ~7 min on Intel(R) Core(TM) i7-8565U CPU @ 1.80GHz
kfinal <- group_assign(
  data = sim_surveys,
  id = "id",
  xcoord = "x",
  ycoord = "y",
  time = names(groupperday),
  group_vector = groupperday,
  method = "hclust")

save(kfinal, file = "Outputs/kfinal1000.RData")
