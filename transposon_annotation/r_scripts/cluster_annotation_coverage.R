rm(list=ls())    #this clears the workspace to make sure no leftover variables are floating around
graphics.off()   #close all graphics windows

library(ggplot2)
library(gridExtra)

cv <- read.table("all_hannuus_220k_samples_by_family_cv.tsv", header=T,sep="\t")
sims_of_220k <- read.table("all_hannuus_220k_samples_by_family.tsv",header=T,sep="\t")

covboxplot <- ggplot(sims_of_220k, aes(x=sims_of_220k$Family, y=sims_of_220k[,2]))
covboxplot + stat_boxplot(geom = "errorbar") + geom_boxplot() + stat_summary(fun.y=mean, geom="line", colour="red", aes(group=1))
covboxplot + scale_y_continuous(name="Genomic Proportion") + opts(axis.text.x = theme_text(angle=90, hjust=1)) + xlab("Family")
covboxplot + xlim("RLG-begi","RLC-amov","L1","Copia","RLC-jiliwu","RLG-rewu","RLG-iketas","Gypsy","RLG-rahi","RLC-X","RLC-suwi","RLG-X",
                  "RLG-wily","RLG-kefe","RLG-teda","RLG-ryse","RLC-ogaow","Kolobok","COP18_I_MT","EnSpm","REP_DE","MuDR")

cvplot <- ggplot(cv, aes(x=seq(1:22),y=cv[,2])) + geom_point() + geom_line() +
cvplot + scale_x_discrete(labels=cv[,1],name="Family") + scale_y_continuous(name="Coefficient of Variation")
cvplot + opts(axis.text.x = theme_text(angle=90, hjust=1))

grid.arrange(cvplot, covboxplot, ncol=1)
