library(picante)

## GLS and PGLS
allsuperfam.dat <- read.table("all_15_species_te-superfamilies_cval_raw_bp_counts_tab.tsv.txt",header=T,sep="\t",row.names=1)
superfam.gls <- gls(GenomeSize ~ Gypsy, data = allsuperfam.dat)
superfam.pgls <- gls(GenomeSize ~ Gypsy, correlation = corBrownian(value = 1, tre), data = allsuperfam.dat)
plot(GenomeSize ~ Gypsy, data = allsuperfam.dat)
abline(coef(superfam.gls), lwd = 2, col = "black")
abline(coef(superfam.pgls), lwd = 2, col = "red")
anova(superfam.pgls)


## Kcalc
tre <- read.tree("ape_tree_15_aster_species")
dat <- read.table("all_15_species_te-families_cval_raw_bp_counts_tab_for_picante.txt",row.names=1,header=T,sep="")
fam_Kcalc_stats <- apply(dat, 2, Kcalc, tre)
fam_multiphylosignal_stats <- multiPhylosignal(dat, multi2di(tre))

superfam.dat <- read.table("all_15_species_te-superfamilies_cval_raw_bp_counts_tab_for_picante.txt",row.names=1,header=T,sep="\t")
superfam_Kcalc_stats <- apply(superfam.dat, 2, Kcalc, tre)
superfam_mulitphylosignal_stats <- multiPhylosignal(superfam.dat, multi2di(tre))