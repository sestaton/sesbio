#library(ggplot2)

Aster_shannon_stats <- read.table("all_15_species_te-families_shannon_diversity_stats_refmt.tsv.txt",header=T,sep="\t")

Aster_shannon_stats.H <- Aster_shannon_stats[Aster_shannon_stats$Diversity_statistic == "H",]
Aster_shannon_stats.J <- Aster_shannon_stats[Aster_shannon_stats$Diversity_statistic == "J",]

p <- ggplot(Aster_shannon_stats.H, aes(x=Species, y=Value, color=Diversity_statistic)) +
     geom_point(color="blue",size=3) +
     theme_bw() +
     coord_flip() +
     scale_x_discrete(limits=c("Nasanthus patagonicus","Fulcaldea stuessyi","Gerbera hybrida",
                               "Carthamus tinctorius","Centrapallus pauciflorus","Taraxacum kok-saghyz",
                               "Senecio vulgaris","Pseudognaphalium helleri","Conoclinium coelestinum",
                               "Phoebanthus tenuifolius","Helianthus porteri","Helianthus verticillatus",
                               "Helianthus niveus ssp. tephrodes","Helianthus annuus","Helianthus argophyllus"),
                      labels=c("Nasanthus patagonicus","Fulcaldea stuessyi","Gerbera hybrida",
                               "Carthamus tinctorius","Centrapalus pauciflorus","Taraxacum kok-saghyz",
                               "Senecio vulgaris","Pseudognaphalium helleri","Conoclinium coelestinum",
                               "Phoebanthus tenuifolius","Helianthus porteri","Helianthus verticillatus",
                               "Helianthus niveus ssp. tephrodes","Helianthus annuus","Helianthus argophyllus")) + ylab("Shannon's diversity (H)") 


p2 <- ggplot(Aster_shannon_stats.J, aes(x=Species, y=Value, color=Diversity_statistic)) +
     geom_point(color="black",size=3) +
     theme_bw() +
     coord_flip() +
     scale_x_discrete(limits=c("Nasanthus patagonicus","Fulcaldea stuessyi","Gerbera hybrida",
                               "Carthamus tinctorius","Centrapallus pauciflorus","Taraxacum kok-saghyz",
                               "Senecio vulgaris","Pseudognaphalium helleri","Conoclinium coelestinum",
                               "Phoebanthus tenuifolius","Helianthus porteri","Helianthus verticillatus",
                               "Helianthus niveus ssp. tephrodes","Helianthus annuus","Helianthus argophyllus"),
                      labels=c("Nasanthus patagonicus","Fulcaldea stuessyi","Gerbera hybrida",
                               "Carthamus tinctorius","Centrapalus pauciflorus","Taraxacum kok-saghyz",
                               "Senecio vulgaris","Pseudognaphalium helleri","Conoclinium coelestinum",
                               "Phoebanthus tenuifolius","Helianthus porteri","Helianthus verticillatus",
                               "Helianthus niveus ssp. tephrodes","Helianthus annuus","Helianthus argophyllus")) + ylab("Shannon's evenness (J)") 

gA <- ggplot_gtable(ggplot_build(p + theme(legend.position="none",
                                           axis.text.x = element_text(size=14),
					   axis.text.y = element_text(size=14),
					   axis.title.y = element_blank(),
                                           axis.title.x = element_text(size=16, color="black"))))

gB <- ggplot_gtable(ggplot_build(p2 + theme(legend.position="none",
					    axis.text.x = element_text(size=14),
					    axis.text.y = element_blank(),
                                            axis.ticks.y = element_blank(),
					    axis.title.y = element_blank(),
                                            axis.title.x = element_text(size=16, color="black"))))

maxWidth = grid::unit.pmax(gA$widths[2:3], gB$widths[2:3])
gA$widths[2:3] <- as.list(maxWidth)
gB$widths[2:3] <- as.list(maxWidth)
print(grid.arrange(gA, gB, ncol=2))

#     #scale_color_manual(name="Diversity_statistic", values = c("H" = "blue", "J" = "black")) +
#     theme(axis.title.y = element_blank(), axis.title.x = element_blank(),
#           axis.text.y = element_text(size=14),
#           axis.text.x = element_text(size=14),
#           axis.title.x = element_text(size=16)) + 
#     guides(color=guide_legend(title="Diversity statistic"))
