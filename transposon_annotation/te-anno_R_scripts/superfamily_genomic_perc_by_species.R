rm(list=ls())    #this clears the workspace to make sure no leftover variables are floating around
graphics.off()   #close all graphics windows

## includes
library(plyr)
library(ggplot2)

## functions
as_percent <- function(x) { 100 * sum(x)}

# Ageratina
ager <- read.table("Ager_1m_corrected_annotation_summary_4-24.tsv",header=T,sep="\t")
ager.dd <- ddply(ager, "Superfamily", summarize, Sum = sum(GPerc_corr), TotPerc = as_percent(GPerc_corr))
Species <- rep("Ager",length(ager.dd[,1]))
ager.dd.df <- cbind(Species, ager.dd)
names(ager.dd.df) <- c("Species","Superfamily","Sum","TotPerc")

# Hann
hann <- read.table("Ann1238_1m_corrected_annotation_summary_4-24.tsv",header=T,sep="\t")
hann.dd <- ddply(hann, "Superfamily", summarize, Sum = sum(GPerc_corr), TotPerc = as_percent(GPerc_corr))
Species <- rep("Hann",length(hann.dd[,1]))
hann.dd.df <- cbind(Species, hann.dd)
names(hann.dd.df) <- c("Species","Superfamily","Sum","TotPerc")

# Gnaph
gnaph <- read.table("Gnaph_1m_corrected_annotation_summary_4-24.tsv",header=T,sep="\t")
gnaph.dd <- ddply(gnaph, "Superfamily", summarize, Sum = sum(GPerc_corr), TotPerc = as_percent(GPerc_corr))
Species <- rep("Gnaph",length(gnaph.dd[,1]))
gnaph.dd.df <- cbind(Species, gnaph.dd)
names(gnaph.dd.df) <- c("Species","Superfamily","Sum","TotPerc")

# Senecio
sene <- read.table("Senecio_1m_corrected_annotation_summary_4-24.tsv",header=T,sep="\t")
sene.dd <- ddply(sene, "Superfamily", summarize, Sum = sum(GPerc_corr), TotPerc = as_percent(GPerc_corr))
Species <- rep("Sene",length(sene.dd[,1]))
sene.dd.df <- cbind(Species, sene.dd)
names(sene.dd.df) <- c("Species","Superfamily","Sum","TotPerc")

#CP
cp <- read.table("CP_1m_corrected_annotation_summary_4-24.tsv",header=T,sep="\t")
cp.dd <- ddply(cp, "Superfamily", summarize, Sum = sum(GPerc_corr), TotPerc = as_percent(GPerc_corr))
Species <- rep("CP",length(cp.dd[,1]))
cp.dd.df <- cbind(Species, cp.dd)
names(cp.dd.df) <- c("Species","Superfamily","Sum","TotPerc")

#TKS
tks <- read.table("TKS_1m_corrected_annotation_summary_4-24.tsv",header=T,sep="\t")
tks.dd <- ddply(tks, "Superfamily", summarize, Sum = sum(GPerc_corr), TotPerc = as_percent(GPerc_corr))
Species <- rep("TKS",length(tks.dd[,1]))
tks.dd.df <- cbind(Species, tks.dd)
names(tks.dd.df) <- c("Species","Superfamily","Sum","TotPerc")

#Saff
saff <- read.table("Saff_1m_corrected_annotation_summary_4-24.tsv",header=T,sep="\t")
saff.dd <- ddply(saff, "Superfamily", summarize, Sum = sum(GPerc_corr), TotPerc = as_percent(GPerc_corr))
Species <- rep("Saff",length(saff.dd[,1]))
saff.dd.df <- cbind(Species, saff.dd)
names(saff.dd.df) <- c("Species","Superfamily","Sum","TotPerc")

#Gerb
gerb <- read.table("Gerb_1m_corrected_annotation_summary_4-24.tsv",header=T,sep="\t")
gerb.dd <- ddply(gerb, "Superfamily", summarize, Sum = sum(GPerc_corr), TotPerc = as_percent(GPerc_corr))
Species <- rep("Gerb",length(gerb.dd[,1]))
gerb.dd.df <- cbind(Species, gerb.dd)
names(gerb.dd.df) <- c("Species","Superfamily","Sum","TotPerc")

#Dasy
dasy <- read.table("Dasy_1m_corrected_annotation_summary_4-24.tsv",header=T,sep="\t")
dasy.dd <- ddply(dasy, "Superfamily", summarize, Sum = sum(GPerc_corr), TotPerc = as_percent(GPerc_corr))
Species <- rep("Dasy",length(dasy.dd[,1]))
dasy.dd.df <- cbind(Species, dasy.dd)
names(dasy.dd.df) <- c("Species","Superfamily","Sum","TotPerc")

# Calyc
calyc <- read.table("Calyc_1m_corrected_annotation_summary_4-24.tsv",header=T,sep="\t")
calyc.dd <- ddply(calyc, "Superfamily", summarize, Sum = sum(GPerc_corr), TotPerc = as_percent(GPerc_corr))
Species <- rep("Caly",length(calyc.dd[,1]))
calyc.dd.df <- cbind(Species, calyc.dd)
names(calyc.dd.df) <- c("Species","Superfamily","Sum","TotPerc")

##merge
all_species_merged.df <- Reduce(function(...) merge(...,all=T), list(ager.dd.df, hann.dd.df, gnaph.dd.df, sene.dd.df,
                                                                     cp.dd.df, tks.dd.df, saff.dd.df, gerb.dd.df, dasy.dd.df,
                                                                     calyc.dd.df))

# copy df & filter by threshold
fullcp <- all_species_merged.df
filtered_full_df <- fullcp[fullcp$TotPerc >= 2.0,]
bp <- ggplot(filtered_full_df, aes(x=Species, y=TotPerc, fill=Superfamily)) + ylim(0,100) + geom_bar(stat="identity",color="black") + coord_flip()
bp + scale_x_discrete(limits=c("Caly","Dasy","Gerb","Saff","TKS","CP","Sene","Gnaph","Hann","Ager"),
                      labels=c("Nasanthus sp.","Dasyphyllum sp.","Gerbera hybrida","Carthamus tinctorius",
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
