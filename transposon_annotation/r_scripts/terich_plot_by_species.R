#library(ggplot2)
#library(plyr)

fam_size_cval_reps <- read.table("Aster_15_species_famct_cval_fammean_all_reps.txt",header=T,sep="\t")
fam_num_sum <- ddply(fam_size_cval_reps, "Species", summarize,
                     FamNum = mean(Families),
                     FamSize = mean(Family_mean),
                     GenomeSize = mean(Genome_size))

### TE richness plot
rich_plot <- ggplot(fam_num_sum, aes(Species, FamNum)) +
  geom_hline(yintercept = mean(fam_num_sum$FamNum), color = "red", linetype = "solid") +
  geom_vline(xintercept = 9.8, color = "black", linetype = "dashed") +
#  geom_rect(data=fam_num_sum, aes(xmin=9.8, xmax=Inf, ymin=-Inf, ymax=mean(fam_num_sum$FamNum)),
#            fill=t, color="black", alpha=0.02, linetype=0) +
  geom_point(size=3) 

### TE family size plot
size_plot <- ggplot(fam_num_sum, aes(Species, FamSize)) +
  geom_hline(yintercept = mean(fam_num_sum$FamSize), color = "red", linetype = "solid") +
  geom_vline(xintercept = 9.8, color = "black", linetype = "dashed") +
#  geom_rect(data=fam_num_sum, aes(xmin=9.8, xmax=Inf, ymin=mean(fam_num_sum$FamSize), ymax=Inf),
#            fill=t, color="black", alpha=0.02, linetype=0) +
  geom_point(size=3)

## align plots
gA <- ggplot_gtable(ggplot_build(rich_plot + scale_x_discrete(limits=c("Calyc","Dasy","Gerb","Saff","CP","TKS","Sene","Gnaph","Ager",
                                      "Phoeb","Hport","Hvert","Hteph","Ann1238","Harg"),
                             labels=c("Nasanthus patagonicus","Fulcaldea stuessyi","Gerbera hybrida",
                                      "Carthamus tinctorius","Centrapalus pauciflorus","Taraxacum kok-saghyz",
                                      "Senecio vulgaris","Pseudognaphalium helleri","Conoclinium coelestinum",
                                      "Phoebanthus tenuifolius","Helianthus porteri","Helianthus verticillatus",
                                      "Helianthus niveus ssp. tephrodes","Helianthus annuus","Helianthus argophyllus")) +
    theme_bw() +
    coord_flip()  +
    ylab("TE family richness") +
    theme(axis.text.x=element_text(size=14, color = "black"),
          axis.text.y=element_text(size=14, color = "black"),
          axis.title.x=element_text(size=16, color = "black"),
          axis.title.y=element_blank(),
          plot.margin = unit(c(1,.5,1,1), "lines"))))

#plot.margin = unit(c(1,0,1,1), "lines")
gB <- ggplot_gtable(ggplot_build(size_plot + scale_x_discrete(limits=c("Calyc","Dasy","Gerb","Saff","CP","TKS","Sene","Gnaph","Ager",
                                      "Phoeb","Hport","Hvert","Hteph","Ann1238","Harg"),
                             labels=c("Nasanthus patagonicus","Fulcaldea stuessyi","Gerbera hybrida",
                                      "Carthamus tinctorius","Centrapalus pauciflorus","Taraxacum kok-saghyz",
                                      "Senecio vulgaris","Pseudognaphalium helleri","Conoclinium coelestinum",
                                      "Phoebanthus tenuifolius","Helianthus porteri","Helianthus verticillatus",
                                      "Helianthus niveus ssp. tephrodes","Helianthus annuus","Helianthus argophyllus")) +
    theme_bw() +
    coord_flip()  +
    ylab("Mean TE family size [%]") +
    theme(axis.text.x=element_text(size=14, color = "black"),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.title.x=element_text(size=16, color = "black"),
          axis.title.y=element_blank(),
          plot.margin = unit(c(1,1,1,.5), "lines"))))

maxWidth = grid::unit.pmax(gA$widths[2:3], gB$widths[2:3])
gA$widths[2:3] <- as.list(maxWidth)
gB$widths[2:3] <- as.list(maxWidth)
print(grid.arrange(gA, gB, nrow=1))
