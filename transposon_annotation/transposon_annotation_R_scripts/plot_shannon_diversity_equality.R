#library(ggplot2)

Aster_shannon_stats <- read.table("all_15_species_te-families_shannon_diversity_stats_refmt.tsv.txt",header=T,sep="\t")

ggplot(Aster_shannon_stats, aes(x=Species, y=Value, color=Diversity_statistic)) +
     geom_point(size=3) +
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
                               "Helianthus niveus ssp. tephrodes","Helianthus annuus","Helianthus argophyllus")) + 
     scale_color_manual(name="Diversity_statistic", values = c("H" = "blue", "J" = "black")) +
     theme(axis.title.y = element_blank(), axis.title.x = element_blank(),
           axis.text.y = element_text(size=14),
           axis.text.x = element_text(size=14),
           axis.title.x = element_text(size=16)) + 
     guides(color=guide_legend(title="Diversity statistic"))