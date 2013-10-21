library(ggplot2)
library(lattice)

# Ager
ager <- read.table("Ageratina_HA383cp.bln",header=F,sep="\t")
Species <- rep("Ageratina",length(ager[,1]))
ager.df <- cbind(Species, ager)
ager_cpbls <- subset(ager.df, select = c(Species, V1, V3))
names(ager_cpbls) <- c("Species","Read","PID")

#Ann
hann <- read.table("Ann1238_HA383cp.bln",header=F,sep="\t")
Species <- rep("H.annuus",length(hann[,1]))
hann.df <- cbind(Species, hann)
hann_cpbls <- subset(hann.df, select = c(Species, V1, V3))
names(hann_cpbls) <- c("Species","Read","PID")

#Dasy
dasy <- read.table("Dasyphyllum_HA383cp.bln",header=F,sep="\t")
Species <- rep("Dasyphyllum",length(dasy[,1]))
dasy.df <- cbind(Species, dasy)
dasy_cpbls <- subset(dasy.df, select = c(Species, V1, V3))
names(dasy_cpbls) <- c("Species","Read","PID")

#gnaph
gnaph <- read.table("Gnaph_HA383cp.bln",header=F,sep="\t")
Species <- rep("Gnaph",length(gnaph[,1]))
gnaph.df <- cbind(Species, gnaph)
gnaph_cpbls <- subset(gnaph.df, select = c(Species, V1, V3))
names(gnaph_cpbls) <- c("Species","Read","PID")

#Hport
hport <- read.table("Hport_HA383cp.bln",header=F,sep="\t")
Species <- rep("H.port.",length(hport[,1]))
hport.df <- cbind(Species, hport)
hport_cpbls <- subset(hport.df, select = c(Species, V1, V3))
names(hport_cpbls) <- c("Species","Read","PID")

#Hvert
hvert <- read.table("Hvert_HA383cp.bln",header=F,sep="\t")
Species <- rep("H.vert.",length(hvert[,1]))
hvert.df <- cbind(Species, hvert)
hvert_cpbls <- subset(hvert.df, select = c(Species, V1, V3))
names(hvert_cpbls) <- c("Species","Read","PID")

#Saff
saff <- read.table("Saff_HA383cp.bln",header=F,sep="\t")
Species <- rep("Carthamus",length(saff[,1]))
saff.df <- cbind(Species, saff)
saff_cpbls <- subset(saff.df, select = c(Species, V1, V3))
names(saff_cpbls) <- c("Species","Read","PID")

#TKS
tks <- read.table("TKS_HA383cp.bln",header=F,sep="\t")
Species <- rep("Taraxacum KS",length(tks[,1]))
tks.df <- cbind(Species, tks)
tks_cpbls <- subset(tks.df, select = c(Species, V1, V3))
names(tks_cpbls) <- c("Species","Read","PID")

#calyc
caly <- read.table("Calyc_HA383cp.bln",header=F,sep="\t")
Species <- rep("Calycera",length(caly[,1]))
caly.df <- cbind(Species, caly)
caly_cpbls <- subset(caly.df, select = c(Species, V1, V3))
names(caly_cpbls) <- c("Species","Read","PID")

#gerb
gerb <- read.table("Gerb_HA383cp.bln",header=F,sep="\t")
Species <- rep("Gerbera",length(gerb[,1]))
gerb.df <- cbind(Species, gerb)
gerb_cpbls <- subset(gerb.df, select = c(Species, V1, V3))
names(gerb_cpbls) <- c("Species","Read","PID")

#Harg
harg <- read.table("Harg_HA383cp.bln",header=F,sep="\t")
Species <- rep("H.arg.",length(harg[,1]))
harg.df <- cbind(Species, harg)
harg_cpbls <- subset(harg.df, select = c(Species, V1, V3))
names(harg_cpbls) <- c("Species","Read","PID")

#Hteph
hteph <- read.table("Hteph_HA383cp.bln",header=F,sep="\t")
Species <- rep("H.teph.",length(hteph[,1]))
hteph.df <- cbind(Species, hteph)
hteph_cpbls <- subset(hteph.df, select = c(Species, V1, V3))
names(hteph_cpbls) <- c("Species","Read","PID")

#Phoeb
phoeb <- read.table("Phoeb_HA383cp.bln",header=F,sep="\t")
Species <- rep("Phoeb",length(phoeb[,1]))
phoeb.df <- cbind(Species, phoeb)
phoeb_cpbls <- subset(phoeb.df, select = c(Species, V1, V3))
names(phoeb_cpbls) <- c("Species","Read","PID")

#sene
sene <- read.table("Senecio_HA383cp.bln",header=F,sep="\t")
Species <- rep("Senecio",length(sene[,1]))
sene.df <- cbind(Species, sene)
sene_cpbls <- subset(sene.df, select = c(Species, V1, V3))
names(sene_cpbls) <- c("Species","Read","PID")

#CP
cp <- read.table("CP_HA383cp.bln",header=F,sep="\t")
Species <- rep("Centrapallus",length(cp[,1]))
cp.df <- cbind(Species, cp)
cp_cpbls <- subset(cp.df, select = c(Species, V1, V3))
names(cp_cpbls) <- c("Species","Read","PID")

all_species_merged.df <- Reduce(function(...) merge(...,all=T), list(ager_cpbls, hann_cpbls, gnaph_cpbls, sene_cpbls,
                                                                     cp_cpbls, tks_cpbls, saff_cpbls, gerb_cpbls, dasy_cpbls,
                                                                     caly_cpbls, harg_cpbls, hport_cpbls, hteph_cpbls, hvert_cpbls, phoeb_cpbls))

pdf("all_species_reads_to_ref_pid.pdf")
xyplot(PID~Species,data=all_species_merged.df, 
    main="Percent identity of all chloroplast reads mapped to HA383cp genome", 
    ylab="Species", xlab="Percent identity")