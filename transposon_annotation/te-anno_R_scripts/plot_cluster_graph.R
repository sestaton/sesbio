library(rgexf)
library(igraph)

source("frlayout.R")

gr <- read.gexf("CL100_graph.gexf")
l <- layout.fruchterman.reingold(g,dim=3,verbose=F)
plot(g, layout=l, vertex.label=NA)