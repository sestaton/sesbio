library(ggplot2)

## plot family-level phylogenetic signal
fam_mulitphylostats <- read.table("all_15_species_fam_multiphylosignal_stats_annot_filtered.txt",header=T,sep="\t")
dp <- ggplot() + geom_point(data=fam_mulitphylostats, aes(Family, K, color=Superfamily)) + theme(axis.text.x=element_text(angle=90, hjust=1, size=6, color="black"))

## plot superfamily-level phylogenetic signal
Kstats_reorder <- read.csv("all_15_species_superfam_multiphylosignal_stats_sig-only.csv")
rectangles <- data.frame( xmin = c(1,5,13), xmax = c(4, 12, 25), ymin = c(0, 0, 0), ymax=c(3.5, 3.5, 3.5), t=c("LTR-RT","non-LTR-RT","Class II")
gp <- ggplot() + geom_point(data=Kstats_reorder, aes(Species, K)) + geom_rect(data=rectangles, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill=t, color="black", alpha=0.2, linetype=0)
                         
gp + scale_x_discrete(limits=c("Copia","ERV1","ERV2","Gypsy",
                               "L1","L2","NeSL","Penelope","R1","R2","RTE",
                               "SINE2.tRNA","CRE","Crack","Daphne","EnSpm",
                               "Harbinger","Helitron","Mariner.Tc1","MuDR",
                               "Polinton","Sola","Tad1","Tx1","hAT"),
                               labels=c("Copia","ERV1","ERV2","Gypsy",
                                        "L1","L2","NeSL","Penelope","R1","R2","RTE",
                                        "SINE2.tRNA","CRE","Crack","Daphne","EnSpm",
                                        "Harbinger","Helitron","Mariner.Tc1","MuDR",
                                        "Polinton","Sola","Tad1","Tx1","hAT"))
