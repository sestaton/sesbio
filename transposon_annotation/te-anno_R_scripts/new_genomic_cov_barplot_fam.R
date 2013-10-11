# This script generates a bar plot for the genomic proportion for a given family
#
# 9/26/13 SES
#
# TODO:

rm(list=ls())    #this clears the workspace to make sure no leftover variables are floating around
graphics.off()   #close all graphics windows

## includes
library(plyr)
library(ggplot2)

## functions
as_percent <- function(x) { 100 * sum(x)}

## read in data
#asters <- read.table("all_10_species_table.tsv",sep="\t",header=T)

## make subsets for each species, and summarize the percent composition
# Ageratina
ager <- read.table("ager_transposome_results_9-10_annotations_summary.tsv", sep="\t", header=T)
ager.summ <- ddply(ager, "Family", summarize, sum = sum(GenomePerc), as_perc = as_percent(GenomePerc), Mbp = sum * 1973)
Species <- rep("Ager",length(ager.summ[,1]))
ager.summ.f <- cbind(Species, ager.summ)
names(ager.summ.f) <- c("Species", "Family","Sum","GenomePerc","Mbp")
rm(ager.summ)

# H. annuus
ann1238 <- read.table("ann1238_transposome_results_9-10_annotations_summary.tsv", sep="\t", header=T)
ann.summ <- ddply(ann1238, "Family", summarize, sum = sum(GenomePerc), as_perc = as_percent(GenomePerc), Mbp = sum * 3600)
Species <- rep("Hann",length(ann.summ[,1]))
ann.summ.f <- cbind(Species, ann.summ)
names(ann.summ.f) <- c("Species", "Family","Sum","GenomePerc","Mbp")
rm(ann.summ)

# H. arg
harg <- read.table("harg_transposome_results_9-10_annotations_summary.tsv", sep="\t", header=T)
harg.summ <- ddply(harg, "Family", summarize, sum = sum(GenomePerc), as_perc = as_percent(GenomePerc), Mbp = sum * 4328)
Species <- rep("Harg",length(harg.summ[,1]))
harg.summ.f <- cbind(Species, harg.summ)
names(harg.summ.f) <- c("Species", "Family","Sum","GenomePerc","Mbp")
rm(harg.summ)

# H. port
hport <- read.table("hport_transposome_results_9-10_annotations_summary.tsv", sep="\t", header=T)
hport.summ <- ddply(hport, "Family", summarize, sum = sum(GenomePerc), as_perc = as_percent(GenomePerc), Mbp = sum * 3400)
Species <- rep("Hport",length(hport.summ[,1]))
hport.summ.f <- cbind(Species, hport.summ)
names(hport.summ.f) <- c("Species", "Family","Sum","GenomePerc","Mbp")
rm(hport.summ)

# H. teph
hteph <- read.table("hteph_transposome_results_9-10_annotations_summary.tsv", sep="\t", header=T)
hteph.summ <- ddply(hteph, "Family", summarize, sum = sum(GenomePerc), as_perc = as_percent(GenomePerc), Mbp = sum * 3570)
Species <- rep("Hteph",length(hteph.summ[,1]))
hteph.summ.f <- cbind(Species, hteph.summ)
names(hteph.summ.f) <- c("Species", "Family","Sum","GenomePerc","Mbp")
rm(hteph.summ)

# H. vert
hvert <- read.table("hvert_transposome_results_9-11_annotations_summary.tsv", sep="\t", header=T)
hvert.summ <- ddply(hvert, "Family", summarize, sum = sum(GenomePerc), as_perc = as_percent(GenomePerc), Mbp = sum * 3500)
Species <- rep("Hvert",length(hvert.summ[,1]))
hvert.summ.f <- cbind(Species, hvert.summ)
names(hvert.summ.f) <- c("Species", "Family","Sum","GenomePerc","Mbp")
rm(hvert.summ)

# Phoeb
phoeb <- read.table("phoeb_transposoem_9-25_cluster_log_annotations_summary.tsv", sep="\t", header=T)
phoeb.summ <- ddply(hvert, "Family", summarize, sum = sum(GenomePerc), as_perc = as_percent(GenomePerc), Mbp = sum * 3500)
Species <- rep("Phoeb",length(phoeb.summ[,1]))
phoeb.summ.f <- cbind(Species, phoeb.summ)
names(phoeb.summ.f) <- c("Species", "Family","Sum","GenomePerc","Mbp")
rm(phoeb.summ)

# CP
cp <- read.table("cp_transposome_results_9-10_annotations_summary.tsv", sep="\t", header=T)
cp.summ <- ddply(cp, "Family", summarize, sum = sum(GenomePerc), as_perc = as_percent(GenomePerc), Mbp = sum * 1702)
Species <- rep("CP",length(cp.summ[,1]))
cp.summ.f <- cbind(Species, cp.summ)
names(cp.summ.f) <- c("Species", "Family","Sum","GenomePerc","Mbp")
rm(cp.summ)

# Calyc
calyc <- read.table("calyc_transposome_results_1m_annotations_summary.tsv", sep="\t", header=T)
calyc.summ <- ddply(calyc, "Family", summarize, sum = sum(GenomePerc), as_perc = as_percent(GenomePerc), Mbp = sum * 1887)
Species <- rep("Calyc",length(calyc.summ[,1]))
calyc.summ.f <- cbind(Species, calyc.summ)
names(calyc.summ.f) <- c("Species", "Family","Sum","GenomePerc","Mbp")
rm(calyc.summ)

# Gerb
gerb <- read.table("gerb_transposome_results_9-10_annotations_summary.tsv", sep="\t", header=T)
gerb.summ <- ddply(gerb, "Family", summarize, sum = sum(GenomePerc), as_perc = as_percent(GenomePerc), Mbp = sum * 2494)
Species <- rep("Gerb",length(gerb.summ[,1]))
gerb.summ.f <- cbind(Species, gerb.summ)
names(gerb.summ.f) <- c("Species", "Family","Sum","GenomePerc","Mbp")
rm(gerb.summ)

# gnaph
gnaph <- read.table("gnaph_transposome_results_9-10_annotations_summary.tsv", sep="\t", header=T)
gnaph.summ <- ddply(gnaph, "Family", summarize, sum = sum(GenomePerc), as_perc = as_percent(GenomePerc), Mbp = sum * 1249)
Species <- rep("Gnaph",length(gnaph.summ[,1]))
gnaph.summ.f <- cbind(Species, gnaph.summ)
names(gnaph.summ.f) <- c("Species", "Family","Sum","GenomePerc","Mbp")
rm(gnaph.summ)

# Saff
saff <- read.table("saff_transposome_results_9-10_annotations_summary.tsv", sep="\t", header=T)
saff.summ <- ddply(saff, "Family", summarize, sum = sum(GenomePerc), as_perc = as_percent(GenomePerc), Mbp = sum * 1364)
Species <- rep("Saff",length(saff.summ[,1]))
saff.summ.f <- cbind(Species, saff.summ)
names(saff.summ.f) <- c("Species", "Family","Sum","GenomePerc","Mbp")
rm(saff.summ)

# Senecio
sene <- read.table("sene_transposome_results_annotations_summary.tsv", sep="\t", header=T)
sene.summ <- ddply(sene, "Family", summarize, sum = sum(GenomePerc), as_perc = as_percent(GenomePerc), Mbp = sum * 1540)
Species <- rep("Sene",length(sene.summ[,1]))
sene.summ.f <- cbind(Species, sene.summ)
names(sene.summ.f) <- c("Species", "Family","Sum","GenomePerc","Mbp")
rm(sene.summ)

# TKS
tks <- read.table("tks_transposome_results_annotations_summary.tsv", sep="\t", header=T)
tks.summ <- ddply(tks, "Family", summarize, sum = sum(GenomePerc), as_perc = as_percent(GenomePerc), Mbp = sum * 1871)
Species <- rep("TKS",length(tks.summ[,1]))
tks.summ.f <- cbind(Species, tks.summ)
names(tks.summ.f) <- c("Species", "Family","Sum","GenomePerc","Mbp")
rm(tks.summ)

## merge data
all_species_merged.df <- Reduce(function(...) merge(...,all=T), list(ager.summ.f, phoeb.summ.f, ann.summ.f, harg.summ.f, hport.summ.f, hteph.summ.f,
                                                                     hvert.summ.f, gnaph.summ.f, sene.summ.f, cp.summ.f, tks.summ.f,
                                                                     saff.summ.f, gerb.summ.f, calyc.summ.f))

# filter the data
fullcp <- all_species_merged.df
filtered_full_df <- fullcp[fullcp$GenomePerc >= 1.9,]

## plot

bp <- ggplot(filtered_full_df, aes(x=Species, y=Mbp, fill=Family)) + geom_bar(stat="identity",color="black") + coord_flip() + scale_x_discrete(limits=c("Calyc","Gerb","Saff","TKS","CP","Sene","Gnaph","Hvert","Hteph","Hport","Harg","Hann","Phoeb","Ager"),
                      labels=c("Nasanthus patagonicus",
                        "Gerbera hybrida",
                        "Carthamus tinctorius",
                        "Taraxacum kok-saghyz",
                        "Centrapallus pauciflorus",
                        "Senecio vulgaris",
                        "Pseudognaphalium helleri",
                        "Helianthus arghophyllus",
                        "Helianthus porteri",
                        "Helianthus tephrodes",
                        "Helianthus verticillatus",
                        "Helianthus annuus",
                        "Phoebathus sp.",
                        "Ageratina altissima"))

#+ scale_fill_manual(values = c("Gypsy" = "darkgreen",
#                      "Copia" = "aquamarine4",
#                      "DIRS" = "lightgreen",
#                      "L1" = "darkgrey",
#                      "EnSpm" = "azure2",
#                      "Crypton" = "chartreuse",
#                      "hAT" = "darkkhaki",
#                      "Helitron" = "darkolivegreen3",
#                      "Mariner/Tc1" = "chocolate",
#                      "L2" = "aquamarine",
#                      "MuDR" = "dimgrey",
#                      "ERV1" = "darkblue",
#                      "CR1" = "chartreuse3",
#                      "Penelope" = "cyan",
#                      "R1" = "black")) + xlab(label = "Mbp") + theme(
#                                                axis.text.x = element_text(color = "black"),                                                               #                                                 axis.text.y = element_text(color = "black", size = 14),
#                                                axis.title.x = element_text(color = "black", size = 14),
#                                                axis.title.y = element_blank())
