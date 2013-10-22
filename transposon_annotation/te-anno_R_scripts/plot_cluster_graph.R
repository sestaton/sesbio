# 1. Take an edge list from a cluster and convert it to GEXF format
# 2. Load the frlayout.R script from fgclust for plotting
# 3. Parse the graph with the rgexf library
# 4. Plot the graph with igraph

library(rgexf)
library(igraph)

source("frlayout.R")

gr <- read.gexf("CL100_graph.gexf")
l <- layout.fruchterman.reingold(g,dim=3,verbose=F)
plot(g, layout=l, vertex.label=NA)