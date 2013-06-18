# This script generates a line plot over a boxplot showing the coefficient
# of variation (top) a boxplot of genomic proportion for a given family
#
# 6/17/13 SES
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
ager <- read.table("Ageratina_500k_merged_cluster_rep_6-17_annotations_summary_edit.tsv", sep="\t", header=T)
ager.summ <- ddply(ager, "Superfamily", summarize, sum = sum(GPerc), as_perc = as_percent(GPerc), Mbp = sum * 1973)
Species <- rep("Ager",length(ager.summ[,1]))
ager.summ.f <- cbind(Species, ager.summ)
names(ager.summ.f) <- c("Species", "Superfamily","Sum","GPerc","Mbp")
rm(ager.summ)

# H. annuus
ann1238 <- read.table("Ann1238_500k_merged_cluster_rep_6-17_annotations_summary_edit.tsv", sep="\t", header=T)
ann.summ <- ddply(ann1238, "Superfamily", summarize, sum = sum(GPerc), as_perc = as_percent(GPerc), Mbp = sum * 3600)
Species <- rep("Ann",length(ann.summ[,1]))
ann.summ.f <- cbind(Species, ann.summ)
names(ann.summ.f) <- c("Species", "Superfamily","Sum","GPerc","Mbp")
rm(ann.summ)

# CP
cp <- read.table("CP_500k_merged_cluster_rep_6-17_annotations_summary_edit.tsv", sep="\t", header=T)
cp.summ <- ddply(cp, "Superfamily", summarize, sum = sum(GPerc), as_perc = as_percent(GPerc), Mbp = sum * 1702)
Species <- rep("CP",length(cp.summ[,1]))
cp.summ.f <- cbind(Species, cp.summ)
names(cp.summ.f) <- c("Species", "Superfamily","Sum","GPerc","Mbp")
rm(cp.summ)

# Calyc
calyc <- read.table("Calyc_500k_merged_cluster_rep_6-17_annotations_summary_edit.tsv", sep="\t", header=T)
calyc.summ <- ddply(calyc, "Superfamily", summarize, sum = sum(GPerc), as_perc = as_percent(GPerc), Mbp = sum * 1887)
Species <- rep("Calyc",length(calyc.summ[,1]))
calyc.summ.f <- cbind(Species, calyc.summ)
names(calyc.summ.f) <- c("Species", "Superfamily","Sum","GPerc","Mbp")
rm(calyc.summ)

# Dasy
dasy <- read.table("Dasy_500k_merged_cluster_rep_6-17_annotations_summary_edit.tsv", sep="\t", header=T)
dasy.summ <- ddply(dasy, "Superfamily", summarize, sum = sum(GPerc), as_perc = as_percent(GPerc), Mbp = sum * 2533)
Species <- rep("Dasy",length(dasy.summ[,1]))
dasy.summ.f <- cbind(Species, dasy.summ)
names(dasy.summ.f) <- c("Species", "Superfamily","Sum","GPerc","Mbp")
rm(dasy.summ)

# Gerb
gerb <- read.table("Gerb_500k_merged_cluster_rep_6-17_annotations_summary.tsv", sep="\t", header=T)
gerb.summ <- ddply(gerb, "Superfamily", summarize, sum = sum(GPerc), as_perc = as_percent(GPerc), Mbp = sum * 2494)
Species <- rep("Gerb",length(gerb.summ[,1]))
gerb.summ.f <- cbind(Species, gerb.summ)
names(gerb.summ.f) <- c("Species", "Superfamily","Sum","GPerc","Mbp")
rm(gerb.summ)

# gnaph
gnaph <- read.table("Gnaph_500k_merged_cluster_rep_6-17_annotations_summary_edit.tsv", sep="\t", header=T)
gnaph.summ <- ddply(gnaph, "Superfamily", summarize, sum = sum(GPerc), as_perc = as_percent(GPerc), Mbp = sum * 1249)
Species <- rep("Gnaph",length(gnaph.summ[,1]))
gnaph.summ.f <- cbind(Species, gnaph.summ)
names(gnaph.summ.f) <- c("Species", "Superfamily","Sum","GPerc","Mbp")
rm(gnaph.summ)

# Saff
saff <- read.table("Saff_500k_merged_cluster_rep_6-17_annotations_summary_edit.tsv", sep="\t", header=T)
saff.summ <- ddply(saff, "Superfamily", summarize, sum = sum(GPerc), as_perc = as_percent(GPerc), Mbp = sum * 1364)
Species <- rep("Saff",length(saff.summ[,1]))
saff.summ.f <- cbind(Species, saff.summ)
names(saff.summ.f) <- c("Species", "Superfamily","Sum","GPerc","Mbp")
rm(saff.summ)

# Senecio
sene <- read.table("Senecio_500k_merged_cluster_rep_6-17_annotations_summary_edit.tsv", sep="\t", header=T)
sene.summ <- ddply(sene, "Superfamily", summarize, sum = sum(GPerc), as_perc = as_percent(GPerc), Mbp = sum * 1540)
Species <- rep("Sene",length(sene.summ[,1]))
sene.summ.f <- cbind(Species, sene.summ)
names(sene.summ.f) <- c("Species", "Superfamily","Sum","GPerc","Mbp")
rm(sene.summ)

# TKS
tks <- read.table("TKS_500k_merged_cluster_rep_6-17_annotations_summary_edit.tsv", sep="\t", header=T)
tks.summ <- ddply(tks, "Superfamily", summarize, sum = sum(GPerc), as_perc = as_percent(GPerc), Mbp = sum * 1871)
Species <- rep("TKS",length(tks.summ[,1]))
tks.summ.f <- cbind(Species, tks.summ)
names(tks.summ.f) <- c("Species", "Superfamily","Sum","GPerc","Mbp")
rm(tks.summ)

## merge data
all_species_merged.df <- Reduce(function(...) merge(...,all=T), list(ager.summ.f, ann.summ.f, gnaph.summ.f, sene.summ.f,
                                                                     cp.summ.f, tks.summ.f, saff.summ.f, gerb.summ.f, dasy.summ.f,
                                                                     calyc.summ.f))

# filter the data
fullcp <- all_species_merged.df
filtered_full_df <- fullcp[fullcp$GPerc >= 0.2,]

## plot
#bp <- ggplot(filtered_full_df, aes(x=Species, y=Mbp, fill=Superfamily)) + geom_bar(stat="identity")  + coord_flip() + scale_x_discrete(limits=c("Calyc","Dasy","Gerb","Saff","TKS","CP","Sene","Gnaph","Ann","Ager"),
#                      labels=c("Nasanthus sp.","Dasyphyllum sp.","Gerbera hybrida","Carthamus tinctorius",
#                               "Taraxacum kok-saghyz","Centrapallus pauciflorus","Senecio vulgaris","Gnaph. sp.",
#                               "Helianthus annuus","Ageratina sp."))
#+ scale_y_discrete(limits=c("0","500","1500","2000","2500"),labels=c("0","1500","2000","2500","3000","3500"))

bp <- ggplot(filtered_full_df, aes(x=Species, y=Mbp, fill=Superfamily)) + geom_bar(stat="identity",color="black") + coord_flip() + scale_x_discrete(limits=c("Calyc","Dasy","Gerb","Saff","TKS","CP","Sene","Gnaph","Ann","Ager"),
                      labels=c("Nasanthus patagonicus","Fulcaldea stuessyi","Gerbera hybrida","Carthamus tinctorius",
                               "Taraxacum kok-saghyz","Centrapallus pauciflorus","Senecio vulgaris","Gnaph. sp.",
                               "Helianthus annuus","Ageratina sp.")) + scale_fill_manual(values = c("Copia" = "aquamarine4", "DIRS" = "lightgreen",
                                                                       "EnSpm" = "azure2", "Gypsy" = "darkgreen",
                                                                       "Harbinger" = "chartreuse", "hAT" = "darkkhaki",
                                                                       "Helitron" = "darkolivegreen3","Kiri" = "chocolate",
                                                                        "L2" = "aquamarine","MuDR" = "dimgrey",
                                                                        "SINE" = "darkblue","CR1" = "chartreuse3",
                                                                        "Polinton" = "cyan4","R1" = "black","L1" = "cadetblue3",
                                                                        "Crack" = "darkgrey")) + theme(axis.text.x = element_text(color = "black"),
                                                                                             axis.text.y = element_text(color = "black", size = 16),
                                                                                             axis.title.x = element_blank(), axis.title.y = element_blank())
