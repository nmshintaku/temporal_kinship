#Supplemental - SEAMAMMS Visuals

library(SocGen)
library(igraph)
library(dplyr)
library(ggplot2)

combined_pairs <- read.csv("Outputs/combined_pairs.csv")
sightings <- read.csv("Outputs/sightings.csv")
dates <- sort(unique(sightings$Observation.Date))

##########
#Boxplot##
##########

ggplot(combined_pairs, aes(x = interaction(State1, State2), y = weight)) +
  geom_boxplot() +
  labs(x = "Repro States", y = "Association Coefficient") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


##########
#Networks#
##########

#Create network from combined pairs df

network <- graph_from_data_frame(combined_pairs[,c("ID1", "ID2", "weight", "DyadML", "Closekin")], 
                                 directed = FALSE)

network <- set_vertex_attr(network, "ReproCat",
                           value = sightings$ReproCat[match(V(network)$name, sightings$Repro.ID)])

synchrony_graph <- subgraph_from_edges(network, eids = E(network)[E(network)$Closekin == "Y"])


#Plot network for all females
V(network)$color <- "skyblue"
plot(network, vertex.label = NA, vertex.size = 5, edge.arrow.size = 0.3, vertex.color = V(network)$color)

#Plot network for females who are both pregnant
preg_ver <- V(network)[V(network)$ReproCat == "preg"]
preg_subgraph <- induced_subgraph(network, preg_ver)
plot(preg_subgraph, vertex.label = NA, vertex.size = 10, edge.arrow.size = 0.5, vertex.color = "deeppink4")

#Plot network for females who are both lactating 
lact_ver <- V(network)[V(network)$ReproCat == "lact"]
lact_subgraph <- induced_subgraph(network, lact_ver)
plot(lact_subgraph, vertex.label = NA, vertex.size = 10, edge.arrow.size = 0.5, vertex.color = "coral2")


#Plot network for all females with close kin edge weighted
#E(synchrony_graph)$width <- 1
#plot(synchrony_graph, vertex.label = NA, vertex.size = 5, edge.arrow.size = 0.3, vertex.color = V(network)$color,
     #edge.width = E(synchrony_graph)$width)


