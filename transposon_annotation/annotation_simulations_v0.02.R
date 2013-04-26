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
library(gridExtra)

## functions
co.var <- function(x) ( 100*sd(x)/mean(x) )
as_percent <- function(x) { 100 * x}

all_sims <- read.table("all_hannuus_full_annotations_summaries_4-21.tsv",sep="\t",header=F)
names(all_sims) <- c("ReadNum","Superfamily","Family","ReadCt/ReadsWithHit","HitPerc","GPerc")

for(i in levels(all_sims$Superfamily)) {
  ## Select one family
  covplot <- paste(i,"_genome_coverage_full",".pdf",sep="")
  pdf(file=covplot) #, height=6.30, width=4.31)
  i <- subset(all_sims, Superfamily == i, select = c(ReadNum, Family, GPerc))
  if (length(i$ReadNum) > 2) {
    i.covar <- ddply(i, "ReadNum", summarize, 
                     sd = sd(GPerc), 
                     mean = mean(GPerc), 
                     covar = co.var(GPerc))

    #ok <- complete.cases(i.covar)
    
    i.covar[is.na(i.covar)] <- 0
    if (length(i.covar[,1]) > 3) {
      ## set up line/box plots
      i.lineplot <- ggplot(i.covar, aes(x=seq(1:length(i.covar[,1])), y=covar)) + geom_line() + geom_point() + scale_x_discrete(labels=i.covar[,1], name="Read Number")  
      
      i.boxplot <- ggplot(i, aes(x=factor(i$ReadNum), y=as_percent(GPerc))) + stat_boxplot(geom = "errorbar") + geom_boxplot() 
      
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
    dev.off() # turn the device off after printing
  }
  else {
    unlink(covplot)
  }
}


