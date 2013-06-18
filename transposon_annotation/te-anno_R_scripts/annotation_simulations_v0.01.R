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

all_sims <- read.table("all_hannuus_new_annotations_summaries.tsv",sep="\t",header=F)
names(all_sims) <- c("ReadNum","Superfamily","Family","ReadCt/AllReads","ReadCt/ReadsWithHit","HitPerc","GPerc")

## Select one family
iketas <- subset(all_sims, Family == "RLG-iketas", select = c(ReadNum, Family, GPerc))
iketas.covar <- ddply(iketas, "ReadNum", summarize, 
                      sd = sd(GPerc), 
                      mean = mean(GPerc), 
                      covar = co.var(GPerc))

## set up line/box plots
iketas.lineplot <- ggplot(iketas.covar, aes(x=seq(1:9), y=covar)) + geom_line() + geom_point() + scale_x_discrete(labels=iketas.covar[,1], name="Read Number")  

iketas.boxplot <- ggplot(iketas, aes(x=factor(iketas$ReadNum), y=as_percent(GPerc))) + stat_boxplot(geom = "errorbar") + geom_boxplot() 

## align plots
gA <- ggplot_gtable(ggplot_build(iketas.lineplot 
	+ ylab("Coefficient of Variation") 
	+ theme(axis.title.x = element_blank(), 
                axis.text.x = element_blank(), 
                axis.ticks.x = element_blank(),
                plot.margin = unit(c(1,1,-0.8,1), "lines")) ))
gB <- ggplot_gtable(ggplot_build(iketas.boxplot 
	+ stat_summary(fun.y=mean, geom="line", colour="red", aes(group=1)) 
	+ xlab("Read Number") 
	+ ylab("Genomic Proportion") 
	+ theme(axis.text.x = element_text(angle=90, hjust=1), 
                plot.margin = unit(c(0.5,1,1,1), "lines"))))
maxWidth = grid::unit.pmax(gA$widths[2:3], gB$widths[2:3])
gA$widths[2:3] <- as.list(maxWidth)
gB$widths[2:3] <- as.list(maxWidth)
grid.arrange(gA, gB, ncol=1)

#hann850 <- read.table("Ann1238_500k_rep_new_annotations_summary.tsv",sep="\t",header=T)
#ddply(hann850, "Superfamily", summarize, sum = sum(HitPerc), as_perc = as_percent(sum(HitPerc)))


