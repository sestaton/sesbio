#library(ggplot2)
#library(plyr)

fam_size_cval_reps <- read.table("Aster_15_species_famct_cval_fammean_all_reps.txt",header=T,sep="\t")
fam_num_sum <- ddply(fam_size_cval_reps, "Species", summarize,
                     FamNum = mean(Families),
                     FamSize = mean(Family_mean),
                     GenomeSize = mean(Genome_size))

##mean(fam_num_sum$FamNum)
#45.06667

rich_plot <- ggplot(fam_num_sum, aes(Species, FamNum)) + geom_hline(yintercept = 45.06667, color = "red", type = 2) + geom_point(size=3) 
rich_plot + scale_x_discrete(limits=c("Calyc","Dasy","Gerb","Saff","CP","TKS","Sene","Gnaph","Ager",
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
          axis.title.y=element_blank())
