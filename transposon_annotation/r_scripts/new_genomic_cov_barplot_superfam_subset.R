# This script generates a histrogram show the % genomic proportion for each superfamily
#
# 12/03/13 SES
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
ager.summ <- ddply(ager, "Superfamily", summarize, sum = sum(GenomePerc), as_perc = as_percent(GenomePerc), Cval = 1746, Mbp = sum * 1746)
Species <- rep("Ager",length(ager.summ[,1]))
ager.summ.f <- cbind(Species, ager.summ)
names(ager.summ.f) <- c("Species", "Superfamily","Sum","GenomePerc","Cval", "Mbp")
rm(ager.summ)

# H. annuus
ann1238 <- read.table("ann1238_transposome_results_9-10_annotations_summary.tsv", sep="\t", header=T)
ann.summ <- ddply(ann1238, "Superfamily", summarize, sum = sum(GenomePerc), as_perc = as_percent(GenomePerc), Cval = 3384, Mbp = sum * 3384)
Species <- rep("Hann",length(ann.summ[,1]))
ann.summ.f <- cbind(Species, ann.summ)
names(ann.summ.f) <- c("Species", "Superfamily","Sum","GenomePerc","Cval","Mbp")
rm(ann.summ)

# CP
cp <- read.table("cp_transposome_results_9-10_annotations_summary.tsv", sep="\t", header=T)
cp.summ <- ddply(cp, "Superfamily", summarize, sum = sum(GenomePerc), as_perc = as_percent(GenomePerc), Cval = 3125, Mbp = sum * 3125)
Species <- rep("CP",length(cp.summ[,1]))
cp.summ.f <- cbind(Species, cp.summ)
names(cp.summ.f) <- c("Species", "Superfamily","Sum","GenomePerc","Cval","Mbp")
rm(cp.summ)

# Calyc
calyc <- read.table("calyc_transposome_results_1m_annotations_summary.tsv", sep="\t", header=T)
calyc.summ <- ddply(calyc, "Superfamily", summarize, sum = sum(GenomePerc), as_perc = as_percent(GenomePerc), Cval = 3962, Mbp = sum * 3962)
Species <- rep("Calyc",length(calyc.summ[,1]))
calyc.summ.f <- cbind(Species, calyc.summ)
names(calyc.summ.f) <- c("Species", "Superfamily","Sum","GenomePerc","Cval","Mbp")
rm(calyc.summ)

# Dasy
dasy <- read.table("dasy_transposome_cluster_report_11-17_annotations_summary.tsv", sep="\t", header=T)
dasy.summ <- ddply(dasy, "Superfamily", summarize, sum = sum(GenomePerc), as_perc = as_percent(GenomePerc), Cval = 4182, Mbp = sum * 4182)
Species <- rep("Dasy",length(dasy.summ[,1]))
dasy.summ.f <- cbind(Species, dasy.summ)
names(dasy.summ.f) <- c("Species", "Superfamily","Sum","GenomePerc","Cval","Mbp")
rm(dasy.summ)

# Gerb
gerb <- read.table("gerb_transposome_results_9-10_annotations_summary.tsv", sep="\t", header=T)
gerb.summ <- ddply(gerb, "Superfamily", summarize, sum = sum(GenomePerc), as_perc = as_percent(GenomePerc), Cval = 3861, Mbp = sum * 3861)
Species <- rep("Gerb",length(gerb.summ[,1]))
gerb.summ.f <- cbind(Species, gerb.summ)
names(gerb.summ.f) <- c("Species", "Superfamily","Sum","GenomePerc","Cval","Mbp")
rm(gerb.summ)

# gnaph
gnaph <- read.table("gnaph_transposome_results_9-10_annotations_summary.tsv", sep="\t", header=T)
gnaph.summ <- ddply(gnaph, "Superfamily", summarize, sum = sum(GenomePerc), as_perc = as_percent(GenomePerc), Cval = 2920, Mbp = sum * 2920)
Species <- rep("Gnaph",length(gnaph.summ[,1]))
gnaph.summ.f <- cbind(Species, gnaph.summ)
names(gnaph.summ.f) <- c("Species", "Superfamily","Sum","GenomePerc","Cval","Mbp")
rm(gnaph.summ)

# Saff
saff <- read.table("saff_transposome_results_9-10_annotations_summary.tsv", sep="\t", header=T)
saff.summ <- ddply(saff, "Superfamily", summarize, sum = sum(GenomePerc), as_perc = as_percent(GenomePerc), Cval = 2405, Mbp = sum * 2405)
Species <- rep("Saff",length(saff.summ[,1]))
saff.summ.f <- cbind(Species, saff.summ)
names(saff.summ.f) <- c("Species", "Superfamily","Sum","GenomePerc","Cval","Mbp")
rm(saff.summ)

# TKS
tks <- read.table("tks_transposome_results_annotations_summary.tsv", sep="\t", header=T)
tks.summ <- ddply(tks, "Superfamily", summarize, sum = sum(GenomePerc), as_perc = as_percent(GenomePerc), Cval = 2582, Mbp = sum * 2582)
Species <- rep("TKS",length(tks.summ[,1]))
tks.summ.f <- cbind(Species, tks.summ)
names(tks.summ.f) <- c("Species", "Superfamily","Sum","GenomePerc","Cval","Mbp")
rm(tks.summ)

## merge data
all_species_merged.df <- Reduce(function(...) merge(...,all=T), list(ager.summ.f, ann.summ.f, gnaph.summ.f, cp.summ.f, tks.summ.f,
                                                                     saff.summ.f, gerb.summ.f, calyc.summ.f, dasy.summ.f))

# filter the data
fullcp <- all_species_merged.df
filtered_full_df <- fullcp[fullcp$GenomePerc >= 0.2,]

## plot

bp <- ggplot(filtered_full_df, aes(x=Species, y=Cval)) + geom_point(aes(size=Cval)) + labs(size="Genome size") + coord_flip() + geom_bar(aes(x=Species, y=Mbp, fill=Superfamily), stat="identity",color="black") + scale_x_discrete(limits=c("Calyc","Dasy", "Gerb","Saff","TKS","CP","Gnaph","Ager", "Hann"), labels=c("Nasanthus patagonicus",
                        "Fulcaldea stuessyi",
                        "Gerbera hybrida",
                        "Carthamus tinctorius",
                        "Taraxacum kok-saghyz",
                        "Centrapallus pauciflorus",
                        "Pseudognaphalium helleri",
                        "Conoclinium coelestinum",
                        "Helianthus annuus")) + scale_fill_manual(name="Superfamily",values = c("Gypsy" = "darkgreen",
                                             "Copia" = "aquamarine4",
                                                                      "DIRS" = "lightgreen",
                                                                      "L1" = "darkgrey",
                                                                      "EnSpm" = "azure2",
                                                                      "Harbinger" = "azure4",
                                                                      "Crypton" = "chartreuse",
                                                                      "hAT" = "darkkhaki",
                                                                      "Helitron" = "darkolivegreen3",
                                                                      "Mariner/Tc1" = "chocolate",
                                                                      "L2" = "aquamarine",
                                                                      "MuDR" = "dimgrey",
                                                                      "ERV1" = "darkblue",
                                                                      "CR1" = "chartreuse3",
                                                                      "Penelope" = "cyan",
                                                                      "R1" = "black")) + theme(axis.text.x = element_text(color = "black"),                                                                                                                                  axis.text.y = element_text(color = "black", size = 14),
                                                                                               axis.title.x = element_blank(),
                                                                                               axis.title.y = element_blank())
