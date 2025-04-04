# Plot an example network from one year

sightings <- read.csv("Outputs/sightings.csv")

sightings <- sightings[which(sightings$ReproSightings >= 10 &
  sightings$ReproCat != "unknown"), ]

sightings$Year <- format(as.Date(sightings$Observation.Date), "%Y")

count <- sightings[!duplicated(sightings[, c("Dolphin.ID", "Year")]), ]

table(count$Year) |> sort()

# select 2004 as the year with the most individuals

sightings2004 <- sightings[which(sightings$Year == 2004), ]

library(SocGen)

network <- simple_ratio(
  sightings = sightings2004,
  group_variable = "Observation.ID",
  dates = "Observation.Date",
  IDs = "Dolphin.ID"
)

states <- sightings2004[
  !duplicated(sightings2004[, c("Dolphin.ID")]),
  c("Dolphin.ID", "ReproCat")
]

colrs <- adjustcolor(hcl.colors(4), alpha.f = 0.8)
basestates <- c("juvenile", "cyc", "preg", "lact")

states$color <- colrs[match(states$ReproCat, basestates)]

library(igraph)

m <- graph_from_adjacency_matrix(network,
  mode = "undirected",
  weighted = TRUE,
  diag = FALSE
)

all(names(V(m)) == states$Dolphin.ID)

V(m)$colr <- states$color

V(m)$degree <- degree(m)

lo <- layout_with_fr(m)

set.seed(1)

pdf(file = "Outputs/example_network.pdf")
plot(m,
  layout = lo,
  vertex.color = V(m)$background,
  vertex.label = NA,
  vertex.size = V(m)$degree / 2 + 3,
  edge.width = edge_attr(m)$weight * 5,
  edge.curved = rep(-.4, length(edge_attr(m)$weight))
)

plot(m,
  add = TRUE,
  layout = lo,
  vertex.color = V(m)$colr,
  vertex.label = NA,
  vertex.size = V(m)$degree / 2 + 3,
  edge.color = NA
)

legend("topright",
  legend = c(
    "juvenile",
    "cycling",
    "pregnant",
    "lactating"
  ),
  pt.bg = colrs, pch = 21, pt.cex = 1.5
)

dev.off()
