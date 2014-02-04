# this script starts with the "filtered_full_df" data frame
# from the script 'new_genomic_cov_barplot_superfam_all_latest.R'

#rm(list=ls())   # in this case, we NEED the left over variable 'filtered_full_df'
graphics.off()   #close all graphics windows

## includes
library(plyr)
library(ggplot2)
library(gridExtra)

filtered_full_df.gypsy <- subset(filtered_full_df, Superfamily == "Gypsy", select = c(Species, Superfamily, GenomePerc))
filtered_full_df.copia <- subset(filtered_full_df, Superfamily == "Copia", select = c(Species, Superfamily, GenomePerc))

## set up dot plots
#color="darkgreen",
p <- ggplot(filtered_full_df.gypsy, aes(x=Species, y=GenomePerc, group=1)) + geom_point(color="darkgreen",size=4) + stat_smooth(method=lm) + scale_x_discrete(limits=c("Calyc","Dasy","Gerb","Saff","CP","TKS","Sene","Gnaph","Ager","Phoeb","Hport","Hvert","Hteph","Hann","Harg"), labels=c("Nasanthus patagonicus","Fulcaldea stuessyi","Gerbera hybrida","Carthamus tinctorius","Centrapalus pauciflorus","Taraxacum kok-saghyz","Senecio vulgaris","Pseudognaphalium helleri","Conoclinium coelestinum","Phoebanthus tenuifolius","Helianthus porteri","Helianthus verticillatus","Helianthus niveus ssp. tephrodes","Helianthus annuus","Helianthus argophyllus"))

#color="aquamarine4",
p2 <- ggplot(filtered_full_df.copia, aes(x=Species, y=GenomePerc, group=1)) + geom_point(color="aquamarine4",size=4) + stat_smooth(method=lm) + scale_x_discrete(limits=c("Calyc","Dasy","Gerb","Saff","CP","TKS","Sene","Gnaph","Ager","Phoeb","Hport","Hvert","Hteph","Hann","Harg"), labels=c("Nasanthus patagonicus","Fulcaldea stuessyi","Gerbera hybrida","Carthamus tinctorius","Centrapalus pauciflorus","Taraxacum kok-saghyz","Senecio vulgaris","Pseudognaphalium helleri","Conoclinium coelestinum","Phoebanthus tenuifolius","Helianthus porteri","Helianthus verticillatus","Helianthus niveus ssp. tephrodes","Helianthus annuus","Helianthus argophyllus"))

## align plots
gA <- ggplot_gtable(ggplot_build(p + coord_flip() + ylab("Genomic Coverage [%]") + theme(axis.text.y = element_text(color = "black", size = 14),
                                                                                         axis.text.x = element_text(color = "black", size = 14),
                                                                                         axis.title.y = element_blank(),
                                                                                         plot.margin = unit(c(1,.5,1,1), "lines"))))
#plot.margin = unit(c(1,0,1,1), "lines")
gB <- ggplot_gtable(ggplot_build(p2 + coord_flip() + ylab("Genomic Coverage [%]") + theme(axis.text.y = element_blank(),
                                                                                          axis.text.x = element_text(color = "black", size = 14),
                                                                                          axis.title.y = element_blank(),
                                                                                          axis.ticks.y = element_blank(),
                                                                                          plot.margin = unit(c(1,1,1,.5), "lines"))))
#plot.margin = unit(c(1,1,1,0), "lines")

maxWidth = grid::unit.pmax(gA$widths[2:3], gB$widths[2:3])
gA$widths[2:3] <- as.list(maxWidth)
gB$widths[2:3] <- as.list(maxWidth)
print(grid.arrange(gA, gB, nrow=1))

