rm(list=ls())    #this clears the workspace to make sure no leftover variables are floating around
graphics.off()   #close all graphics windows

## includes
library(plyr)
library(ggplot2)
library(gridExtra)

## functions
co.var <- function(x) ( 100*sd(x)/mean(x) )

## read data and give columns names
all_sims <- read.table("all_hannuus_new_annotations_summaries.tsv",sep="\t",header=F)
names(all_sims) <- c("ReadNum","Superfamily","Family","ReadCt/AllReads","ReadCt/ReadsWithHit","HitPerc","GPerc")

## select one family
iketas <- subset(all_sims, Family == "RLG-iketas", select = c(ReadNum, Family, HitPerc))
iketas.covar <- ddply(iketas, "ReadNum", summarize, 
                      sd = sd(HitPerc), 
                      mean = mean(HitPerc), 
                      covar = co.var(HitPerc))

## set up plots
iketas.lineplot <- ggplot(iketas.covar, aes(x=seq(1:8), y=covar)) + geom_line() + geom_point()
iketas.lineplot + scale_x_discrete(labels=iketas.covar[,1], name="Read Number") 
iketas.lineplot + ylab("Coefficient of Variation") #scale_y_continuous(name="Coefficient of Variation") + 
iketas.lineplot + theme(axis.title.x = element_blank(), 
                        axis.text.x = element_blank(), 
                        axis.ticks.x = element_blank(),
                        plot.margin = unit(c(1,1,-0.8,1), "lines"))
                                                                                                                                                                                                          plot.margin = unit(c(1,1,-0.8,1), "lines"))
iketas.boxplot <- ggplot(iketas, aes(x=factor(iketas$ReadNum), y=HitPerc)) + stat_boxplot(geom = "errorbar") + geom_boxplot() 
iketas.boxplot + stat_summary(fun.y=mean, geom="line", colour="red", aes(group=1)) 
iketas.boxplot + scale_y_continuous(name="Genomic Proportion") 
iketas.boxplot + theme(axis.text.x = element_text(angle=90, hjust=1), 
                       plot.margin = unit(c(0.1,1,1,1), "lines")) 

## align plots
#gA <- ggplot_gtable(ggplot_build(iketas.lineplot))
#gB <- ggplot_gtable(ggplot_build(iketas.boxplot))
#maxWidth = grid::unit.pmax(gA$widths[2:3], gB$widths[2:3])
#gA$widths[2:3] <- as.list(maxWidth)
#gB$widths[2:3] <- as.list(maxWidth)
#grid.arrange(gA, gB, ncol=1)

grid.arrange(iketas.lineplot + ggtitle("Variation in Genomic Proportion for RLG-iketas"), 
             iketas.boxplot + xlab("Read Number") + ylab("Genomic Proportion") 
             + stat_summary(fun.y=mean, geom="line", colour="red", aes(group=1)) 
             + theme(axis.text.x = element_text(angle=90, hjust=1)), ncol=1)