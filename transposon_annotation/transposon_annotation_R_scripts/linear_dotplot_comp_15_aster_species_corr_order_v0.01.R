# this script starts with the "filtered_full_df" data frame
# from the script 'new_genomic_cov_barplot_superfam_all_latest.R'

#rm(list=ls())   # in this case, we NEED the left over variable 'filtered_full_df'
graphics.off()   #close all graphics windows

## includes
library(plyr)
library(ggplot2)
library(gridExtra)

lm_eqn = function(df){
     m = lm(Gypsy ~ Species, df);
     eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                      list(a = format(coef(m)[1], digits = 2), 
                           b = format(coef(m)[2], digits = 2), 
                           r2 = format(summary(m)$r.squared, digits = 3)))
     as.character(as.expression(eq));                 
}

lm2_eqn = function(df){
     m = lm(Copia ~ Species, df);
     eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,
                      list(a = format(coef(m)[1], digits = 2),
                           b = format(coef(m)[2], digits = 2),
                           r2 = format(summary(m)$r.squared, digits = 3)))
     as.character(as.expression(eq));
}

gyp_cop_perc <- read.table("../all_15_species_gypsy_copia_perc.txt",header=T,sep="\t")

fit1 <- lm(Gypsy ~ Species, data = gyp_cop_perc)
fit2 <- lm(Copia ~ Species, data = gyp_cop_perc)

## set up dot plots
#
p <- ggplot(gyp_cop_perc, aes(x=Species, y=Gypsy, group=1)) +
  geom_point(color="darkgreen",size=4,position="jitter") +
  stat_smooth(method=lm, color="red") + 
  scale_x_discrete(limits=c("Calyc","Dasy","Gerb","Saff",
                            "CP","TKS","Sene","Gnaph","Ager",
                            "Phoeb","Hport","Hvert","Hteph",
                            "Hann","Harg"),
                   labels=c("Nasanthus patagonicus","Fulcaldea stuessyi","Gerbera hybrida",
                            "Carthamus tinctorius","Centrapalus pauciflorus","Taraxacum kok-saghyz",
                            "Senecio vulgaris","Pseudognaphalium helleri","Conoclinium coelestinum",
                            "Phoebanthus tenuifolius","Helianthus porteri","Helianthus verticillatus",
                            "Helianthus niveus ssp. tephrodes","Helianthus annuus","Helianthus argophyllus")) +
  coord_flip() + 
  annotate("text",y=60,x=9,label = lm_eqn(gyp_cop_perc), size = 3, parse=TRUE) +
  theme(axis.text.y = element_text(color = "black", size = 14),
        axis.text.x = element_text(color = "black", size = 14),
        axis.title.y = element_blank())
        #plot.margin = unit(c(1,.5,1,1), "lines"))

#
p2 <- ggplot(gyp_cop_perc, aes(x=Species, y=Copia, group=1)) +
  geom_point(color="aquamarine4",size=4,position="jitter") +
  stat_smooth(method=lm, color="red") + 
  scale_x_discrete(limits=c("Calyc","Dasy","Gerb","Saff",
                            "CP","TKS","Sene","Gnaph","Ager",
                            "Phoeb","Hport","Hvert","Hteph",
                            "Hann","Harg"),
                   labels=c("Nasanthus patagonicus","Fulcaldea stuessyi","Gerbera hybrida",
                            "Carthamus tinctorius","Centrapalus pauciflorus","Taraxacum kok-saghyz",
                            "Senecio vulgaris","Pseudognaphalium helleri","Conoclinium coelestinum",
                            "Phoebanthus tenuifolius","Helianthus porteri","Helianthus verticillatus",
                            "Helianthus niveus ssp. tephrodes","Helianthus annuus","Helianthus argophyllus")) + 
  coord_flip() + 
  annotate("text", y=30, x=9, label = lm2_eqn(gyp_cop_perc), size = 3, parse=TRUE) +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(color = "black", size = 14),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank())
        #plot.margin = unit(c(1,1,1,.5), "lines"))

## align plots
gA <- ggplot_gtable(ggplot_build(p + ylab("Genomic Coverage [%]") + theme(axis.text.y = element_text(color = "black", size = 14),
                                                                          axis.text.x = element_text(color = "black", size = 14),
                                                                          axis.title.y = element_blank(),
                                                                          plot.margin = unit(c(1,.5,1,1), "lines"))))
#plot.margin = unit(c(1,0,1,1), "lines")
gB <- ggplot_gtable(ggplot_build(p2 + ylab("Genomic Coverage [%]") + theme(axis.text.y = element_blank(),
                                                                           axis.text.x = element_text(color = "black", size = 14),
                                                                           axis.title.y = element_blank(),
                                                                           axis.ticks.y = element_blank(),
                                                                           plot.margin = unit(c(1,1,1,.5), "lines"))))

maxWidth = grid::unit.pmax(gA$widths[2:3], gB$widths[2:3])
gA$widths[2:3] <- as.list(maxWidth)
gB$widths[2:3] <- as.list(maxWidth)
print(grid.arrange(gA, gB, nrow=1))
#grid.arrange(p, p2, nrow=1)
