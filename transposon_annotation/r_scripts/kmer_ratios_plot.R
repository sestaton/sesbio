###############################################################
# Plot the uniqueness ratios for a range of k-mers  
# Spencer Evan Staton 
# started 3.13.10
#
# Description: read in a table of k-mer lengths with uniiqueness
# ratios by running: gt tallymer occratio with the 'relative' 
# flag and plot the distribution. 
# 
# TODO: write a script that calls: gt tallymer to get k-mer
# statistics and then writes some results, formatted correctly, that 
# are read and plotted by this R script, like this: 
# R <- read_length_plot.R from the commandline. Will need to write
# the plot to file with pdf()  
#
# 
##############################################################
rm(list=ls()) #this clears the workspace to make sure no leftover variables are floating around
graphics.off(); #close all graphics windows
#library(HH) # this is for writing a png file
library(ggplot2) # this is for making nice plots easier
#require(stats) # for lowess plotting

## Do this manual monkey style df creation for each species
# Ageratina
ager <- read.table("Ageratina_CAGATC_prinseq_trimmed_paired_sub1m_10mer_tallymer_occratio_nohead.txt", sep="", fill=T,header=F)
Species <- rep("Ager",length(ager[,1]))
ager.df <- cbind(Species, ager)
names(ager.df) <- c("Species", "MerLen","Count","Ratio")
rm(ager)

# Nasanthus
calyc <- read.table("Calyc_GATCAG_prinseq_trimmed_paired_sub1m_10mer_tallymer_occratio_nohead.txt", sep="", fill=T,header=F)
Species <- rep("Calyc",length(calyc[,1]))
calyc.df <- cbind(Species, calyc)
names(calyc.df) <- c("Species", "MerLen","Count","Ratio")
rm(calyc)

# Centrapallus
cp <- read.table("CP_AGTCAA_prinseq_trimmed_paired_sub1m_10mer_tallymer_occratio_nohead.txt", sep="", fill=T,header=F)
Species <- rep("CP",length(cp[,1]))
cp.df <- cbind(Species, cp)
names(cp.df) <- c("Species", "MerLen","Count","Ratio")
rm(cp)

# Safflower
saff <- read.table("Saff_GATCAG_prinseq_trimmed_paired_sub1m_10mer_tallymer_occratio_nohead.txt",sep="",fill=T,header=F)
Species <- rep("Saff",length(saff[,1]))
saff.df <- cbind(Species, saff)
names(saff.df) <- c("Species", "MerLen","Count","Ratio")
rm(saff)

# Senecio
sene <- read.table("Senecio_TTAGGC_prinseq_trimmed_paired_sub1m_10mer_tallymer_occratio_nohead.txt",sep="",fill=T,header=F)
Species <- rep("Sene",length(sene[,1]))
sene.df <- cbind(Species, sene)
names(sene.df) <- c("Species", "MerLen","Count","Ratio")
rm(sene)

# Dasy
dasy <- read.table("Dasyphllum_ATCACG_prinseq_trimmed_paired_sub1m_10mer_tallymer_occratio_nohead.txt",sep="",fill=T,header=F)
Species <- rep("Dasy",length(dasy[,1]))
dasy.df <- cbind(Species, dasy)
names(dasy.df) <- c("Species", "MerLen","Count","Ratio")
rm(dasy)

# Gerb
gerb <- read.table("Gerb_GGCTAC_prinseq_trimmed_paired_sub1m_10mer_tallymer_occratio_nohead.txt",sep="",fill=T,header=F)
Species <- rep("Gerb",length(gerb[,1]))
gerb.df <- cbind(Species, gerb)
names(gerb.df) <- c("Species", "MerLen","Count","Ratio")
rm(gerb)

# Gnap
gnap <- read.table("Gnaph_ACAGTG_prinseq_trimmed_paired_sub1m_10mer_tallymer_occratio_nohead.txt",sep="",fill=T,header=F)
Species <- rep("Gnap",length(gnap[,1]))
gnap.df <- cbind(Species, gnap)
names(gnap.df) <- c("Species", "MerLen","Count","Ratio")
rm(gnap)

# Hann
hann <- read.table("PI603989_prinseq_trimmed_paired_sub1m_10mer_tallymer_occratio_nohead.txt",sep="",fill=T,header=F)
Species <- rep("Hann",length(hann[,1]))
hann.df <- cbind(Species, hann)
names(hann.df) <- c("Species", "MerLen","Count","Ratio")
rm(hann)

# TKS
tks <- read.table("TKS_CAGATC_prinseq_trimmed_paired_sub1m_10mer_tallymer_occratio_nohead.txt",sep="",fill=T,header=F)
Species <- rep("TKS",length(tks[,1]))
tks.df <- cbind(Species, tks)
names(tks.df) <- c("Species", "MerLen","Count","Ratio")
rm(tks)

## merge them
ager_sene <- merge(ager.df, sene.df,all=T)
ager_sene_calyc <- merge(ager_sene, calyc.df, all=T);
rm(ager_sene)
ager_sene_calyc_cp <- merge(ager_sene_calyc, cp.df,all=T);
rm(ager_sene_calyc)
ager_sene_calyc_cp_dasy <- merge(ager_sene_calyc_cp, dasy.df,all=T);
rm(ager_sene_calyc_cp)
ager_sene_calyc_cp_dasy_gerb <- merge(ager_sene_calyc_cp_dasy,gerb.df,all=T);
rm(ager_sene_calyc_cp_dasy)
ager_sene_calyc_cp_dasy_gerb_gnaph <- merge(ager_sene_calyc_cp_dasy_gerb,gnap.df,all=T);
rm(ager_sene_calyc_cp_dasy_gerb)
ager_sene_calyc_cp_dasy_gerb_gnaph_hann <- merge(ager_sene_calyc_cp_dasy_gerb_gnaph,hann.df,all=T);
rm(ager_sene_calyc_cp_dasy_gerb_gnaph)
ager_sene_calyc_cp_dasy_gerb_gnaph_hann_sene <- merge(ager_sene_calyc_cp_dasy_gerb_gnaph_hann,sene.df,all=T);
rm(ager_sene_calyc_cp_dasy_gerb_gnaph_hann)
ager_sene_calyc_cp_dasy_gerb_gnaph_hann_sene_tks <- merge(ager_sene_calyc_cp_dasy_gerb_gnaph_hann_sene,tks.df,all=T);
rm(ager_sene_calyc_cp_dasy_gerb_gnaph_hann_sene)

## line plot
ggplot(ager_sene_calyc_cp_dasy_gerb_gnaph_hann_sene_tks, aes(x=MerLen, y=Ratio, color=Species)) + geom_line() + scale_color_brewer(palette="BrBG") + scale_y_continuous(limits=c(0.80,1.00)) + theme(axis.title.x = element_blank(), axis.title.y=element_blank())
