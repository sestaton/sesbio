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
hann_1.2m_s11_sfams <- ddply(hann_1.2m_s11, "Superfamily", summarize, sd = sd(GPerc), mean = mean(GPerc), sum = sum(GPerc))
hann_1.2m_s34_sfams <- ddply(hann_1.2m_s34, "Superfamily", summarize, sd = sd(GPerc), mean = mean(GPerc), sum = sum(GPerc))
hann_1.2m_s56_sfams <- ddply(hann_1.2m_s56, "Superfamily", summarize, sd = sd(GPerc), mean = mean(GPerc), sum = sum(GPerc))
hann_1.2m_merged <- Reduce(function(...) merge(..., all=T), list(hann_1.2m_s11_sfams, hann_1.2m_s34_sfams, hann_1.2m_s56_sfams))
ReadNum <- rep("1220000",length(hann_1.2m_merged[,1]))
hann_1.2m_merged_exp <- cbind(ReadNum, hann_1.2m_merged)

#1.4m
hann_1.4m_s11 <- read.table("hannuus_1.4m_s11_latest_rep_annotations_summary.tsv",sep="\t",header=T)
hann_1.4m_s34 <- read.table("hannuus_1.4m_s34_latest_rep_annotations_summary.tsv",sep="\t",header=T)
hann_1.4m_s56 <- read.table("hannuus_1.4m_s56_latest_rep_annotations_summary.tsv",sep="\t",header=T)
##merge
hann_1.4m_s11_sfams <- ddply(hann_1.4m_s11, "Superfamily", summarize, sd = sd(GPerc), mean = mean(GPerc), sum = sum(GPerc))
hann_1.4m_s34_sfams <- ddply(hann_1.4m_s34, "Superfamily", summarize, sd = sd(GPerc), mean = mean(GPerc), sum = sum(GPerc))
hann_1.4m_s56_sfams <- ddply(hann_1.4m_s56, "Superfamily", summarize, sd = sd(GPerc), mean = mean(GPerc), sum = sum(GPerc))
hann_1.4m_merged <- Reduce(function(...) merge(..., all=T), list(hann_1.4m_s11_sfams, hann_1.4m_s34_sfams, hann_1.4m_s56_sfams))
ReadNum <- rep("1420000",length(hann_1.4m_merged[,1]))
hann_1.4m_merged_exp <- cbind(ReadNum, hann_1.4m_merged)

#1.6m
hann_1.6m_s11 <- read.table("hannuus_1.6m_s11_latest_rep_annotations_summary.tsv",sep="\t",header=T)
hann_1.6m_s34 <- read.table("hannuus_1.6m_s34_latest_rep_annotations_summary.tsv",sep="\t",header=T)
hann_1.6m_s56 <- read.table("hannuus_1.6m_s56_latest_rep_annotations_summary.tsv",sep="\t",header=T)
##merge
hann_1.6m_s11_sfams <- ddply(hann_1.6m_s11, "Superfamily", summarize, sd = sd(GPerc), mean = mean(GPerc), sum = sum(GPerc))
hann_1.6m_s34_sfams <- ddply(hann_1.6m_s34, "Superfamily", summarize, sd = sd(GPerc), mean = mean(GPerc), sum = sum(GPerc))
hann_1.6m_s56_sfams <- ddply(hann_1.6m_s56, "Superfamily", summarize, sd = sd(GPerc), mean = mean(GPerc), sum = sum(GPerc))
hann_1.6m_merged <- Reduce(function(...) merge(..., all=T), list(hann_1.6m_s11_sfams, hann_1.6m_s34_sfams, hann_1.6m_s56_sfams))
ReadNum <- rep("1620000",length(hann_1.6m_merged[,1]))
hann_1.6m_merged_exp <- cbind(ReadNum, hann_1.6m_merged)

#1.8m
hann_1.8m_s11 <- read.table("hannuus_1.8m_s11_latest_rep_annotations_summary.tsv",sep="\t",header=T)
hann_1.8m_s34 <- read.table("hannuus_1.8m_s34_latest_rep_annotations_summary.tsv",sep="\t",header=T)
hann_1.8m_s56 <- read.table("hannuus_1.8m_s56_latest_rep_annotations_summary.tsv",sep="\t",header=T)
##merge
hann_1.8m_s11_sfams <- ddply(hann_1.8m_s11, "Superfamily", summarize, sd = sd(GPerc), mean = mean(GPerc), sum = sum(GPerc))
hann_1.8m_s34_sfams <- ddply(hann_1.8m_s34, "Superfamily", summarize, sd = sd(GPerc), mean = mean(GPerc), sum = sum(GPerc))
hann_1.8m_s56_sfams <- ddply(hann_1.8m_s56, "Superfamily", summarize, sd = sd(GPerc), mean = mean(GPerc), sum = sum(GPerc))
hann_1.8m_merged <- Reduce(function(...) merge(..., all=T), list(hann_1.8m_s11_sfams, hann_1.8m_s34_sfams, hann_1.8m_s56_sfams))
ReadNum <- rep("1820000",length(hann_1.8m_merged[,1]))
hann_1.8m_merged_exp <- cbind(ReadNum, hann_1.8m_merged)

#1m
hann_1m_s11 <- read.table("hannuus_1m_s11_latest_rep_annotations_summary.tsv",sep="\t",header=T)
hann_1m_s34 <- read.table("hannuus_1m_s34_latest_rep_annotations_summary.tsv",sep="\t",header=T)
hann_1m_s56 <- read.table("hannuus_1m_s56_latest_rep_annotations_summary.tsv",sep="\t",header=T)
##merge
hann_1m_s11_sfams <- ddply(hann_1m_s11, "Superfamily", summarize, sd = sd(GPerc), mean = mean(GPerc), sum = sum(GPerc))
hann_1m_s34_sfams <- ddply(hann_1m_s34, "Superfamily", summarize, sd = sd(GPerc), mean = mean(GPerc), sum = sum(GPerc))
hann_1m_s56_sfams <- ddply(hann_1m_s56, "Superfamily", summarize, sd = sd(GPerc), mean = mean(GPerc), sum = sum(GPerc))
hann_1m_merged <- Reduce(function(...) merge(..., all=T), list(hann_1m_s11_sfams, hann_1m_s34_sfams, hann_1m_s56_sfams))
ReadNum <- rep("1020000",length(hann_1m_merged[,1]))
hann_1m_merged_exp <- cbind(ReadNum, hann_1m_merged)

#20k
hann_20k_s11 <- read.table("hannuus_20k_s11_latest_rep_annotations_summary.tsv",sep="\t",header=T)
hann_20k_s34 <- read.table("hannuus_20k_s34_latest_rep_annotations_summary.tsv",sep="\t",header=T)
hann_20k_s56 <- read.table("hannuus_20k_s56_latest_rep_annotations_summary.tsv",sep="\t",header=T)
##merge
hann_20k_s11_sfams <- ddply(hann_20k_s11, "Superfamily", summarize, sd = sd(GPerc), mean = mean(GPerc), sum = sum(GPerc))
hann_20k_s34_sfams <- ddply(hann_20k_s34, "Superfamily", summarize, sd = sd(GPerc), mean = mean(GPerc), sum = sum(GPerc))
hann_20k_s56_sfams <- ddply(hann_20k_s56, "Superfamily", summarize, sd = sd(GPerc), mean = mean(GPerc), sum = sum(GPerc))
hann_20k_merged <- Reduce(function(...) merge(..., all=T), list(hann_20k_s11_sfams, hann_20k_s34_sfams, hann_20k_s56_sfams))
ReadNum <- rep("20000",length(hann_20k_merged[,1]))
hann_20k_merged_exp <- cbind(ReadNum, hann_20k_merged)

#220k
hann_220k_s11 <- read.table("hannuus_220k_s11_latest_rep_annotations_summary.tsv",sep="\t",header=T)
hann_220k_s34 <- read.table("hannuus_220k_s34_latest_rep_annotations_summary.tsv",sep="\t",header=T)
hann_220k_s56 <- read.table("hannuus_220k_s56_latest_rep_annotations_summary.tsv",sep="\t",header=T)
##merge
hann_220k_s11_sfams <- ddply(hann_220k_s11, "Superfamily", summarize, sd = sd(GPerc), mean = mean(GPerc), sum = sum(GPerc))
hann_220k_s34_sfams <- ddply(hann_220k_s34, "Superfamily", summarize, sd = sd(GPerc), mean = mean(GPerc), sum = sum(GPerc))
hann_220k_s56_sfams <- ddply(hann_220k_s56, "Superfamily", summarize, sd = sd(GPerc), mean = mean(GPerc), sum = sum(GPerc))
hann_220k_merged <- Reduce(function(...) merge(..., all=T), list(hann_220k_s11_sfams, hann_220k_s34_sfams, hann_220k_s56_sfams))
ReadNum <- rep("220000",length(hann_220k_merged[,1]))
hann_220k_merged_exp <- cbind(ReadNum, hann_220k_merged)

#420k
hann_420k_s11 <- read.table("hannuus_420k_s11_latest_rep_annotations_summary.tsv",sep="\t",header=T)
hann_420k_s34 <- read.table("hannuus_420k_s34_latest_rep_annotations_summary.tsv",sep="\t",header=T)
hann_420k_s56 <- read.table("hannuus_420k_s56_latest_rep_annotations_summary.tsv",sep="\t",header=T)
##merge
hann_420k_s11_sfams <- ddply(hann_420k_s11, "Superfamily", summarize, sd = sd(GPerc), mean = mean(GPerc), sum = sum(GPerc))
hann_420k_s34_sfams <- ddply(hann_420k_s34, "Superfamily", summarize, sd = sd(GPerc), mean = mean(GPerc), sum = sum(GPerc))
hann_420k_s56_sfams <- ddply(hann_420k_s56, "Superfamily", summarize, sd = sd(GPerc), mean = mean(GPerc), sum = sum(GPerc))
hann_420k_merged <- Reduce(function(...) merge(..., all=T), list(hann_420k_s11_sfams, hann_420k_s34_sfams, hann_420k_s56_sfams))
ReadNum <- rep("420000",length(hann_420k_merged[,1]))
hann_420k_merged_exp <- cbind(ReadNum, hann_420k_merged)

#620k
hann_620k_s11 <- read.table("hannuus_620k_s11_latest_rep_annotations_summary.tsv",sep="\t",header=T)
hann_620k_s34 <- read.table("hannuus_620k_s34_latest_rep_annotations_summary.tsv",sep="\t",header=T)
hann_620k_s56 <- read.table("hannuus_620k_s56_latest_rep_annotations_summary.tsv",sep="\t",header=T)
##merge
hann_620k_s11_sfams <- ddply(hann_620k_s11, "Superfamily", summarize, sd = sd(GPerc), mean = mean(GPerc), sum = sum(GPerc))
hann_620k_s34_sfams <- ddply(hann_620k_s34, "Superfamily", summarize, sd = sd(GPerc), mean = mean(GPerc), sum = sum(GPerc))
hann_620k_s56_sfams <- ddply(hann_620k_s56, "Superfamily", summarize, sd = sd(GPerc), mean = mean(GPerc), sum = sum(GPerc))
hann_620k_merged <- Reduce(function(...) merge(..., all=T), list(hann_620k_s11_sfams, hann_620k_s34_sfams, hann_620k_s56_sfams))
ReadNum <- rep("620000",length(hann_620k_merged[,1]))
hann_620k_merged_exp <- cbind(ReadNum, hann_620k_merged)

#820k
hann_820k_s11 <- read.table("hannuus_820k_s11_latest_rep_annotations_summary.tsv",sep="\t",header=T)
hann_820k_s34 <- read.table("hannuus_820k_s34_latest_rep_annotations_summary.tsv",sep="\t",header=T)
hann_820k_s56 <-read.table("hannuus_820k_s56_latest_rep_annotations_summary.tsv",sep="\t",header=T)
## merge
hann_820k_s11_sfams <- ddply(hann_820k_s11, "Superfamily", summarize, sd = sd(GPerc), mean = mean(GPerc), sum = sum(GPerc))
hann_820k_s34_sfams <- ddply(hann_820k_s34, "Superfamily", summarize, sd = sd(GPerc), mean = mean(GPerc), sum = sum(GPerc))
hann_820k_s56_sfams <- ddply(hann_820k_s56, "Superfamily", summarize, sd = sd(GPerc), mean = mean(GPerc), sum = sum(GPerc))
hann_820k_merged <- Reduce(function(...) merge(..., all=T), list(hann_820k_s11_sfams, hann_820k_s34_sfams, hann_820k_s56_sfams))
ReadNum <- rep("820000",length(hann_820k_merged[,1]))
hann_820k_merged_exp <- cbind(ReadNum, hann_820k_merged)

#
# merge all
#
all_merged <- Reduce(function(...) merge(..., all=T), list(hann_20k_merged_exp, hann_220k_merged_exp, hann_420k_merged_exp, hann_620k_merged_exp,
                                                hann_820k_merged_exp, hann_1m_merged_exp, hann_1.2m_merged_exp, hann_1.4m_merged_exp, hann_1.6m_merged_exp))

#
# plot all
#
for (i in levels(all_merged$Superfamily)) {
  covplot <- paste(i,"_superfamily_genome_coverage_full",".pdf",sep="")
  pdf(file=covplot) #, height=6.30, width=4.31)
  i.sf <- subset(all_merged, Superfamily == i, select = c(ReadNum, sum))
  #p <- ggplot(i.sf, aes(x=ReadNum, y=sum)) + stat_boxplot(geom="errorbar") + geom_boxplot()
  #print(p)
  i.covar <- ddply(i.sf, "ReadNum", summarize, 
                   sd = sd(sum), 
                   mean = mean(sum), 
                   covar = co.var(sum))

  i.covar[is.na(i.covar)] <- 0
  if (length(i.covar[,1]) > 1) {
    i.lineplot <- ggplot(i.covar, aes(x=seq(1:length(i.covar[,1])), y=covar)) + geom_line() + geom_point() + scale_x_discrete(labels=i.covar[,1], name="Read Number")  
      
    i.boxplot <- ggplot(i.sf, aes(x=factor(i.sf$ReadNum), y=sum*100)) + stat_boxplot(geom = "errorbar") + geom_boxplot() 
  
    ## align plots
    gA <- ggplot_gtable(ggplot_build(i.lineplot 
                                     + ylab("Coefficient of Variation") 
                                     + theme(axis.title.x = element_blank(), 
                                             axis.text.x = element_blank(), 
                                             axis.ticks.x = element_blank(),
                                             plot.margin = unit(c(1,1,-0.8,1), "lines")) ))
    gB <- ggplot_gtable(ggplot_build(i.boxplot 
                                     + stat_summary(fun.y=mean, geom="line", colour="red", aes(group=1)) 
                                     + xlab("Read Number") 
                                     + ylab("Genomic Proportion") 
                                     + theme(axis.text.x = element_text(angle=90, hjust=1), 
                                             plot.margin = unit(c(0.5,1,1,1), "lines"))))
    maxWidth = grid::unit.pmax(gA$widths[2:3], gB$widths[2:3])
    gA$widths[2:3] <- as.list(maxWidth)
    gB$widths[2:3] <- as.list(maxWidth)
    print(grid.arrange(gA, gB, ncol=1))
  }
  dev.off()
}
