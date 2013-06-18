# this script starts with the "filtered_full_df" data frame
# from the script 'new_genomic_cov_barplot_superfam_v0.01.R'

#rm(list=ls())   # in this case, we NEED the left over variable 'filtered_full_df'
graphics.off()   #close all graphics windows

## includes
library(plyr)
library(ggplot2)
library(gridExtra)

filtered_full_df.gypsy <- subset(filtered_full_df, Superfamily == "Gypsy", select = c(Species, Superfamily, GPerc))
filtered_full_df.copia <- subset(filtered_full_df, Superfamily == "Copia", select = c(Species, Superfamily, GPerc))

## set up dot plots
p <- ggplot(filtered_full_df.gypsy, aes(x=Species, y=GPerc,group=1)) + stat_smooth(method=lm) + geom_point(color="darkgreen",size=3) + coord_flip() + scale_x_discrete(limits=c("Calyc","Dasy","Gerb","Saff","TKS","CP","Sene","Gnaph","Ann","Ager"),
           labels=c("Nasanthus patagonicus","Fulcaldea stuessyi",
             "Gerbera hybrida","Carthamus tinctorius",
             "Taraxacum kok-saghyz","Centrapallus pauciflorus",
             "Senecio vulgaris","Gnaph. sp.",
             "Helianthus annuus","Ageratina sp."))


p2 <- ggplot(filtered_full_df.copia, aes(x=Species, y=GPerc,group=1)) + stat_smooth(method=lm) + geom_point(color="aquamarine4",size=3) + coord_flip() + scale_x_discrete(limits=c("Calyc","Dasy","Gerb","Saff","TKS","CP","Sene","Gnaph","Ann","Ager"),
              labels=c("Nasanthus patagonicus","Fulcaldea stuessyi",
                "Gerbera hybrida","Carthamus tinctorius",
                "Taraxacum kok-saghyz","Centrapallus pauciflorus",
                "Senecio vulgaris","Gnaph. sp.",
                "Helianthus annuus","Ageratina sp."))

## align plots
gA <- ggplot_gtable(ggplot_build(p 
                                 + ylab("Genomic Proportion") 
                                 #+ theme(axis.title.x = element_blank(), 
                                 #        axis.text.x = element_blank(), 
                                 #        axis.ticks.x = element_blank(),
                                 #        plot.margin = unit(c(1,1,-0.8,1), "lines")) )
))
gB <- ggplot_gtable(ggplot_build(p2 
                                 + ylab("Genomic Proportion") + theme(axis.text.y = element_blank(),
                                                                      axis.title.y = element_blank(),
                                                                      axis.ticks.y = element_blank())))
#                                         plot.margin = unit(c(0.5,1,1,1), "lines"))))


maxWidth = grid::unit.pmax(gA$widths[2:3], gB$widths[2:3])
gA$widths[2:3] <- as.list(maxWidth)
gB$widths[2:3] <- as.list(maxWidth)
print(grid.arrange(gA, gB, nrow=1))

