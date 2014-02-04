## commands for producing a nice looking plot
## NB: the collapse.singles function may not be necessary, depending on the tree

library(ape)
aster.tr <- read.tree("asteraceae_rooted_phylo_nobrlen_latest.newick")
plot.phylo(collapse.singles(aster.tr), edge.width=2,label.offset=0.2,cex=0.5,no.margin=T)