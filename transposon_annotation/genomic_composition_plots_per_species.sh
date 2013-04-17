# This script generates a line plot over a boxplot showing the coefficient
# of variation (top) a boxplot of genomic proportion for a given family
#
# 4/17/13 SES
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
asters <- read.table("all_10_species_table.tsv",sep="\t",header=T)

## make subsets for each species, and summarize the percent composition
# Ageratina
ager <- subset(asters, Species == "Ager", select = c(Species, Superfamily, Family, HitPerc))
ager.summ <- ddply(ager, "Superfamily", summarize, sum = sum(HitPerc), as_perc = as_percent(HitPerc))
Species <- rep("Ager",length(ager.summ[,1]))
ager.summ.f <- cbind(Species, ager.summ)
names(ager.summ.f) <- c("Species", "Superfamily","Sum","GPerc")
rm(ager.summ)

# Calycera (Nasanthus sp.)
calyc <- subset(asters, Species == "Calyc", select = c(Species, Superfamily, Family, HitPerc))
calyc.summ <- ddply(calyc, "Superfamily", summarize, sum = sum(HitPerc), as_perc = as_percent(HitPerc))
Species <- rep("Calyc",length(calyc.summ[,1]))
calyc.summ.f <- cbind(Species, calyc.summ)
names(calyc.summ.f) <- c("Species", "Superfamily","Sum","GPerc")
rm(calyc.summ)

# Centrapallus
cp <- subset(asters, Species == "CP", select = c(Species, Superfamily, Family, HitPerc))
cp.summ <- ddply(cp, "Superfamily", summarize, sum = sum(HitPerc), as_perc = as_percent(HitPerc))
Species <- rep("CP",length(cp.summ[,1]))
cp.summ.f <- cbind(Species, cp.summ)
names(cp.summ.f) <- c("Species", "Superfamily","Sum","GPerc")
rm(cp.summ)

# Dasyphyllum
dasy <- subset(asters, Species == "Dasy", select = c(Species, Superfamily, Family, HitPerc))
dasy.summ <- ddply(dasy, "Superfamily", summarize, sum = sum(HitPerc), as_perc = as_percent(HitPerc))
Species <- rep("Dasy",length(dasy.summ[,1]))
dasy.summ.f <- cbind(Species, dasy.summ)
names(dasy.summ.f) <- c("Species", "Superfamily","Sum","GPerc")
rm(dasy.summ)

# Gerbera hybrida
gerb <- subset(asters, Species == "Gerb", select = c(Species, Superfamily, Family, HitPerc))
gerb.summ <- ddply(gerb, "Superfamily", summarize, sum = sum(HitPerc), as_perc = as_percent(HitPerc))
Species <- rep("Gerb",length(gerb.summ[,1]))
gerb.summ.f <- cbind(Species, gerb.summ)
names(gerb.summ.f) <- c("Species", "Superfamily","Sum","GPerc")
rm(gerb.summ)

# Gnaph sp.
gnap <- subset(asters, Species == "Gnaph", select = c(Species, Superfamily, Family, HitPerc))
gnap.summ <- ddply(gnap, "Superfamily", summarize, sum = sum(HitPerc), as_perc = as_percent(HitPerc))
Species <- rep("Gnaph",length(gnap.summ[,1]))
gnap.summ.f <- cbind(Species, gnap.summ)
names(gnap.summ.f) <- c("Species", "Superfamily","Sum","GPerc")
rm(gnap.summ)

# Helianthus annuus
hann <- subset(asters, Species == "Hann", select = c(Species, Superfamily, Family, HitPerc))
hann.summ <- ddply(hann, "Superfamily", summarize, sum = sum(HitPerc), as_perc = as_percent(HitPerc))
Species <- rep("Hann",length(hann.summ[,1]))
hann.summ.f <- cbind(Species, hann.summ)
names(hann.summ.f) <- c("Species", "Superfamily","Sum","GPerc")
rm(hann.summ)

# Cartamus tinctorius
saff <- subset(asters, Species == "Saff", select = c(Species, Superfamily, Family, HitPerc))
saff.summ <- ddply(saff, "Superfamily", summarize, sum = sum(HitPerc), as_perc = as_percent(HitPerc))
Species <- rep("Saff",length(saff.summ[,1]))
saff.summ.f <- cbind(Species, saff.summ)
names(saff.summ.f) <- c("Species", "Superfamily","Sum","GPerc")
rm(saff.summ)

# Senecio vulgaris
sene <- subset(asters, Species == "Sene", select = c(Species, Superfamily, Family, HitPerc))
sene.summ <- ddply(sene, "Superfamily", summarize, sum = sum(HitPerc), as_perc = as_percent(HitPerc))
Species <- rep("Sene",length(sene.summ[,1]))
sene.summ.f <- cbind(Species, sene.summ)
names(sene.summ.f) <- c("Species", "Superfamily","Sum","GPerc")
rm(sene.summ)

## Taraxacum k-s
tks <- subset(asters, Species == "TKS", select = c(Species, Superfamily, Family, HitPerc))
tks.summ <- ddply(tks, "Superfamily", summarize, sum = sum(HitPerc), as_perc = as_percent(HitPerc))
Species <- rep("TKS",length(tks.summ[,1]))
tks.summ.f <- cbind(Species, tks.summ)
names(tks.summ.f) <- c("Species", "Superfamily","Sum","GPerc")
rm(tks.summ)

## merge data
ager_sene <- merge(ager.summ.f, sene.summ.f,all=T)
ager_sene_calyc <- merge(ager_sene, calyc.summ.f, all=T); rm(ager_sene)
ager_sene_calyc_cp <- merge(ager_sene_calyc, cp.summ.f,all=T); rm(ager_sene_calyc)
ager_sene_calyc_cp_dasy <- merge(ager_sene_calyc_cp, dasy.summ.f,all=T); rm(ager_sene_calyc_cp)
ager_sene_calyc_cp_dasy_gerb <- merge(ager_sene_calyc_cp_dasy,gerb.summ.f,all=T); rm(ager_sene_calyc_cp_dasy)
ager_sene_calyc_cp_dasy_gerb_gnaph <- merge(ager_sene_calyc_cp_dasy_gerb,gnap.summ.f,all=T); rm(ager_sene_calyc_cp_dasy_gerb)
ager_sene_calyc_cp_dasy_gerb_gnaph_hann <- merge(ager_sene_calyc_cp_dasy_gerb_gnaph,hann.summ.f,all=T); rm(ager_sene_calyc_cp_dasy_gerb_gnaph)
ager_sene_calyc_cp_dasy_gerb_gnaph_hann_sene <- merge(ager_sene_calyc_cp_dasy_gerb_gnaph_hann,sene.summ.f,all=T); rm(ager_sene_calyc_cp_dasy_gerb_gnaph_hann)
ager_sene_calyc_cp_dasy_gerb_gnaph_hann_sene_tks <- merge(ager_sene_calyc_cp_dasy_gerb_gnaph_hann_sene,tks.summ.f,all=T); rm(ager_sene_calyc_cp_dasy_gerb_gnaph_hann_sene)

# filter the data
fullcp <- ager_sene_calyc_cp_dasy_gerb_gnaph_hann_sene_tks
filtered_full_df <- fullcp[fullcp$GPerc >= 0.2,]

## plot
bp <- ggplot(all_species_filtered, aes(x=Species, y=GPerc, fill=Superfamily)) + geom_bar(stat="identity")  + coord_flip()
bp + scale_x_discrete(limits=c("Calyc","Dasy","Gerb","Saff","TKS","CP","Sene","Gnaph","Hann","Ager"),
                      labels=c("Nasanthus sp.","Dasyphyllum sp.","Gerbera hybrida","Carthamus tinctorius",
			       "Taraxacum kok-saghyz","Centrapallus pauciflorus","Senecio vulgaris","Gnaph. sp.",
			       "Helianthus annuus","Ageratina sp.")) 
