rm(list=ls())
graphics.off()

library(ggplot2)
library(gridExtra)
library(plyr)


## functions
co.var <- function(x) ( 100*sd(x)/mean(x) )
#as_percent <- function(x) { 100 * x}

#1.2m
hann_1.2m_s11 <- read.table("hannuus_1.2m_s11_latest_rep_annotations_summary.tsv",sep="\t",header=T)
hann_1.2m_s34 <- read.table("hannuus_1.2m_s34_latest_rep_annotations_summary.tsv",sep="\t",header=T)
hann_1.2m_s56 <- read.table("hannuus_1.2m_s56_latest_rep_annotations_summary.tsv",sep="\t",header=T)
##merge
hann_1.2m_merged <- Reduce(function(...) merge(..., all=T), list(hann_1.2m_s11, hann_1.2m_s34, hann_1.2m_s56))
hann_1.2m_fams <- ddply(hann_1.2m_merged, "Family", summarize, sd = sd(GPerc), mean = mean(GPerc), co.var = co.var(GPerc))

ReadNum <- rep("1220000",length(hann_1.2m_fams[,1]))
hann_1.2m_fams_exp <- cbind(ReadNum, hann_1.2m_fams)
hann_1.2m_fams_comp <- hann_1.2m_fams_exp[complete.cases(hann_1.2m_fams_exp),]

#1.4m
hann_1.4m_s11 <- read.table("hannuus_1.4m_s11_latest_rep_annotations_summary.tsv",sep="\t",header=T)
hann_1.4m_s34 <- read.table("hannuus_1.4m_s34_latest_rep_annotations_summary.tsv",sep="\t",header=T)
hann_1.4m_s56 <- read.table("hannuus_1.4m_s56_latest_rep_annotations_summary.tsv",sep="\t",header=T)
##merge
hann_1.4m_merged <- Reduce(function(...) merge(..., all=T), list(hann_1.4m_s11, hann_1.4m_s34, hann_1.4m_s56))
hann_1.4m_fams <- ddply(hann_1.4m_merged, "Family", summarize, sd = sd(GPerc), mean = mean(GPerc), co.var = co.var(GPerc))

ReadNum <- rep("1420000",length(hann_1.4m_fams[,1]))
hann_1.4m_fams_exp <- cbind(ReadNum, hann_1.4m_fams)
hann_1.4m_fams_comp <- hann_1.4m_fams_exp[complete.cases(hann_1.4m_fams_exp),]

#1.6m
hann_1.6m_s11 <- read.table("hannuus_1.6m_s11_latest_rep_annotations_summary.tsv",sep="\t",header=T)
hann_1.6m_s34 <- read.table("hannuus_1.6m_s34_latest_rep_annotations_summary.tsv",sep="\t",header=T)
hann_1.6m_s56 <- read.table("hannuus_1.6m_s56_latest_rep_annotations_summary.tsv",sep="\t",header=T)
##merge
hann_1.6m_merged <- Reduce(function(...) merge(..., all=T), list(hann_1.6m_s11, hann_1.6m_s34, hann_1.6m_s56))
hann_1.6m_fams <- ddply(hann_1.6m_merged, "Family", summarize, sd = sd(GPerc), mean = mean(GPerc), co.var = co.var(GPerc))

ReadNum <- rep("1620000",length(hann_1.6m_fams[,1]))
hann_1.6m_fams_exp <- cbind(ReadNum, hann_1.6m_fams)
hann_1.6m_fams_comp <- hann_1.6m_fams_exp[complete.cases(hann_1.6m_fams_exp),]

#1.8m
hann_1.8m_s11 <- read.table("hannuus_1.8m_s11_latest_rep_annotations_summary.tsv",sep="\t",header=T)
hann_1.8m_s34 <- read.table("hannuus_1.8m_s34_latest_rep_annotations_summary.tsv",sep="\t",header=T)
hann_1.8m_s56 <- read.table("hannuus_1.8m_s56_latest_rep_annotations_summary.tsv",sep="\t",header=T)
##merge
hann_1.8m_merged <- Reduce(function(...) merge(..., all=T), list(hann_1.8m_s11, hann_1.8m_s34, hann_1.8m_s56))
hann_1.8m_fams <- ddply(hann_1.8m_merged, "Family", summarize, sd = sd(GPerc), mean = mean(GPerc), co.var = co.var(GPerc))

ReadNum <- rep("1820000",length(hann_1.8m_fams[,1]))
hann_1.8m_fams_exp <- cbind(ReadNum, hann_1.8m_fams)
hann_1.8m_fams_comp <- hann_1.8m_fams_exp[complete.cases(hann_1.8m_fams_exp),]

#1m
hann_1m_s11 <- read.table("hannuus_1m_s11_latest_rep_annotations_summary.tsv",sep="\t",header=T)
hann_1m_s34 <- read.table("hannuus_1m_s34_latest_rep_annotations_summary.tsv",sep="\t",header=T)
hann_1m_s56 <- read.table("hannuus_1m_s56_latest_rep_annotations_summary.tsv",sep="\t",header=T)
##merge
hann_1m_merged <- Reduce(function(...) merge(..., all=T), list(hann_1m_s11, hann_1m_s34, hann_1m_s56))
hann_1m_fams <- ddply(hann_1m_merged, "Family", summarize, sd = sd(GPerc), mean = mean(GPerc), co.var = co.var(GPerc))

ReadNum <- rep("1000000",length(hann_1m_fams[,1]))
hann_1m_fams_exp <- cbind(ReadNum, hann_1m_fams)
hann_1m_fams_comp <- hann_1m_fams_exp[complete.cases(hann_1m_fams_exp),]

#20k
hann_20k_s11 <- read.table("hannuus_20k_s11_latest_rep_annotations_summary.tsv",sep="\t",header=T)
hann_20k_s34 <- read.table("hannuus_20k_s34_latest_rep_annotations_summary.tsv",sep="\t",header=T)
hann_20k_s56 <- read.table("hannuus_20k_s56_latest_rep_annotations_summary.tsv",sep="\t",header=T)
##merge
hann_20k_merged <- Reduce(function(...) merge(..., all=T), list(hann_20k_s11, hann_20k_s34, hann_20k_s56))
hann_20k_fams <- ddply(hann_20k_merged, "Family", summarize, sd = sd(GPerc), mean = mean(GPerc), co.var = co.var(GPerc))

ReadNum <- rep("20000",length(hann_20k_fams[,1]))
hann_20k_fams_exp <- cbind(ReadNum, hann_20k_fams)
hann_20k_fams_comp <- hann_20k_fams_exp[complete.cases(hann_20k_fams_exp),]

#220k
hann_220k_s11 <- read.table("hannuus_220k_s11_latest_rep_annotations_summary.tsv",sep="\t",header=T)
hann_220k_s34 <- read.table("hannuus_220k_s34_latest_rep_annotations_summary.tsv",sep="\t",header=T)
hann_220k_s56 <- read.table("hannuus_220k_s56_latest_rep_annotations_summary.tsv",sep="\t",header=T)
##merge
hann_220k_merged <- Reduce(function(...) merge(..., all=T), list(hann_220k_s11, hann_220k_s34, hann_220k_s56))
hann_220k_fams <- ddply(hann_220k_merged, "Family", summarize, sd = sd(GPerc), mean = mean(GPerc), co.var = co.var(GPerc))

ReadNum <- rep("220000",length(hann_220k_fams[,1]))
hann_220k_fams_exp <- cbind(ReadNum, hann_220k_fams)
hann_220k_fams_comp <- hann_220k_fams_exp[complete.cases(hann_220k_fams_exp),]

#420k
hann_420k_s11 <- read.table("hannuus_420k_s11_latest_rep_annotations_summary.tsv",sep="\t",header=T)
hann_420k_s34 <- read.table("hannuus_420k_s34_latest_rep_annotations_summary.tsv",sep="\t",header=T)
hann_420k_s56 <- read.table("hannuus_420k_s56_latest_rep_annotations_summary.tsv",sep="\t",header=T)
##merge
hann_420k_merged <- Reduce(function(...) merge(..., all=T), list(hann_420k_s11, hann_420k_s34, hann_420k_s56))
hann_420k_fams <- ddply(hann_420k_merged, "Family", summarize, sd = sd(GPerc), mean = mean(GPerc), co.var = co.var(GPerc))

ReadNum <- rep("420000",length(hann_420k_fams[,1]))
hann_420k_fams_exp <- cbind(ReadNum, hann_420k_fams)
hann_420k_fams_comp <- hann_420k_fams_exp[complete.cases(hann_420k_fams_exp),]

#620k
hann_620k_s11 <- read.table("hannuus_620k_s11_latest_rep_annotations_summary.tsv",sep="\t",header=T)
hann_620k_s34 <- read.table("hannuus_620k_s34_latest_rep_annotations_summary.tsv",sep="\t",header=T)
hann_620k_s56 <- read.table("hannuus_620k_s56_latest_rep_annotations_summary.tsv",sep="\t",header=T)
##merge
hann_620k_merged <- Reduce(function(...) merge(..., all=T), list(hann_620k_s11, hann_620k_s34, hann_620k_s56))
hann_620k_fams <- ddply(hann_620k_merged, "Family", summarize, sd = sd(GPerc), mean = mean(GPerc), co.var = co.var(GPerc))

ReadNum <- rep("620000",length(hann_620k_fams[,1]))
hann_620k_fams_exp <- cbind(ReadNum, hann_620k_fams)
hann_620k_fams_comp <- hann_620k_fams_exp[complete.cases(hann_620k_fams_exp),]

#820k
hann_820k_s11 <- read.table("hannuus_820k_s11_latest_rep_annotations_summary.tsv",sep="\t",header=T)
hann_820k_s34 <- read.table("hannuus_820k_s34_latest_rep_annotations_summary.tsv",sep="\t",header=T)
hann_820k_s56 <- read.table("hannuus_820k_s56_latest_rep_annotations_summary.tsv",sep="\t",header=T)
##merge
hann_820k_merged <- Reduce(function(...) merge(..., all=T), list(hann_820k_s11, hann_820k_s34, hann_820k_s56))
hann_820k_fams <- ddply(hann_820k_merged, "Family", summarize, sd = sd(GPerc), mean = mean(GPerc), co.var = co.var(GPerc))

ReadNum <- rep("820000",length(hann_820k_fams[,1]))
hann_820k_fams_exp <- cbind(ReadNum, hann_820k_fams)
hann_820k_fams_comp <- hann_820k_fams_exp[complete.cases(hann_820k_fams_exp),]

#
# merge all
#
all_merged <- Reduce(function(...) merge(..., all=T), list(hann_20k_fams_comp, hann_220k_fams_comp, hann_420k_fams_comp, hann_620k_fams_comp,
                                                hann_820k_fams_comp, hann_1m_fams_comp, hann_1.2m_fams_comp, hann_1.4m_fams_comp, hann_1.6m_fams_comp))

#
# plot all
#
ggplot(all_merged, aes(mean*100, co.var)) +
  geom_point(aes(alpha=ReadNum)) +
  xlim(0, 5) +
  ylim(-7, 60) +
  stat_smooth(method="loess",fill="blue",alpha=0.1) +
  geom_hline(yintercept = 15, color="green") +
  geom_vline(xintercept = 1, color="green") +
  theme_bw()
