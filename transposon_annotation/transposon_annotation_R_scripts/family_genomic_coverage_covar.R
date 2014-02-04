rm(list=ls())    #this clears the workspace to make sure no leftover variables are floating around
graphics.off()   #close all graphics windows

## includes
library(plyr)
library(ggplot2)
library(gridExtra)

## functions
co.var <- function(x) ( 100*sd(x)/mean(x) )
as_percent <- function(x) { 100 * x}


all_20k <- read.table("all_20k.tsv",header=F,sep="\t")
names(all_20k) <- c("ReadNum","Superfamily", "Family","ReadswithHit", "hitperc", "GPerc")
all_20k.dd <- ddply(all_20k, "Family", summarize, covar = co.var(hitperc), mean = mean(hitperc*100))
all_20k.dd.nna <- na.omit(all_20k.dd)
ReadNum <- rep("20000",length(all_20k.dd.nna[,1]))
all_20k.dd.nna.df <- cbind(ReadNum,all_20k.dd.nna)

all_220k <- read.table("all_220k.tsv",header=F,sep="\t")
names(all_220k) <- c("ReadNum","Superfamily", "Family","ReadswithHit", "hitperc", "GPerc")
all_220k.dd <- ddply(all_220k, "Family", summarize, covar = co.var(hitperc), mean = mean(hitperc*100))
all_220k.dd.nna <- na.omit(all_220k.dd)
ReadNum <- rep("220000",length(all_220k.dd.nna[,1]))
all_220k.dd.nna.df <- cbind(ReadNum,all_220k.dd.nna)


all_420k <- read.table("all_420k.tsv",header=F,sep="\t")
names(all_420k) <- c("ReadNum","Superfamily", "Family","ReadswithHit", "hitperc", "GPerc")
all_420k.dd <- ddply(all_420k, "Family", summarize, covar = co.var(hitperc), mean = mean(hitperc*100))
all_420k.dd.nna <- na.omit(all_420k.dd)
ReadNum <- rep("420000",length(all_420k.dd.nna[,1]))
all_420k.dd.nna.df <- cbind(ReadNum,all_420k.dd.nna)

all_620k <- read.table("all_620k.tsv",header=F,sep="\t")
names(all_620k) <- c("ReadNum","Superfamily", "Family","ReadswithHit", "hitperc", "GPerc")
all_620k.dd <- ddply(all_620k, "Family", summarize, covar = co.var(hitperc), mean = mean(hitperc*100))
all_620k.dd.nna <- na.omit(all_620k.dd)
ReadNum <- rep("620000",length(all_620k.dd.nna[,1]))
all_620k.dd.nna.df <- cbind(ReadNum,all_620k.dd.nna)

all_820k <- read.table("all_820k.tsv",header=F,sep="\t")
names(all_820k) <- c("ReadNum","Superfamily", "Family","ReadswithHit", "hitperc", "GPerc")
all_820k.dd <- ddply(all_820k, "Family", summarize, covar = co.var(hitperc), mean = mean(hitperc*100))
all_820k.dd.nna <- na.omit(all_820k.dd)
ReadNum <- rep("820000",length(all_820k.dd.nna[,1]))
all_820k.dd.nna.df <- cbind(ReadNum,all_820k.dd.nna)

all_1m <- read.table("all_1m.tsv",header=F,sep="\t")
names(all_1m) <- c("ReadNum","Superfamily", "Family","ReadswithHit", "hitperc", "GPerc")
all_1m.dd <- ddply(all_1m, "Family", summarize, covar = co.var(hitperc), mean = mean(hitperc*100))
all_1m.dd.nna <- na.omit(all_1m.dd)
ReadNum <- rep("1020000",length(all_1m.dd.nna[,1]))
all_1m.dd.nna.df <- cbind(ReadNum,all_1m.dd.nna)

all_1.2m <- read.table("all_1.2m.tsv",header=F,sep="\t")
names(all_1.2m) <- c("ReadNum","Superfamily", "Family","ReadswithHit", "hitperc", "GPerc")
all_1.2m.dd <- ddply(all_1.2m, "Family", summarize, covar = co.var(hitperc), mean = mean(hitperc*100))
all_1.2m.dd.nna <- na.omit(all_1.2m.dd)
ReadNum <- rep("1220000",length(all_1.2m.dd.nna[,1]))
all_1.2m.dd.nna.df <- cbind(ReadNum,all_1.2m.dd.nna)

all_1.4m <- read.table("all_1.4m.tsv",header=F,sep="\t")
names(all_1.4m) <- c("ReadNum","Superfamily", "Family","ReadswithHit", "hitperc", "GPerc")
all_1.4m.dd <- ddply(all_1.4m, "Family", summarize, covar = co.var(hitperc), mean = mean(hitperc*100))
all_1.4m.dd.nna <- na.omit(all_1.4m.dd)
ReadNum <- rep("1420000",length(all_1.4m.dd.nna[,1]))
all_1.4m.dd.nna.df <- cbind(ReadNum,all_1.4m.dd.nna)

all_1.6m <- read.table("all_1.6m.tsv",header=F,sep="\t")
names(all_1.6m) <- c("ReadNum","Superfamily", "Family","ReadswithHit", "hitperc", "GPerc")
all_1.6m.dd <- ddply(all_1.6m, "Family", summarize, covar = co.var(hitperc), mean = mean(hitperc*100))
all_1.6m.dd.nna <- na.omit(all_1.6m.dd)
ReadNum <- rep("1620000",length(all_1.6m.dd.nna[,1]))
all_1.6m.dd.nna.df <- cbind(ReadNum,all_1.6m.dd.nna)

all_1.8m <- read.table("all_1.8m.tsv",header=F,sep="\t")
names(all_1.8m) <- c("ReadNum","Superfamily", "Family","ReadswithHit", "hitperc", "GPerc")
all_1.8m.dd <- ddply(all_1.8m, "Family", summarize, covar = co.var(hitperc), mean = mean(hitperc*100))
all_1.8m.dd.nna <- na.omit(all_1.8m.dd)
ReadNum <- rep("1820000",length(all_1.8m.dd.nna[,1]))
all_1.8m.dd.nna.df <- cbind(ReadNum,all_1.8m.dd.nna)
all_merged_finally <- Reduce(function(...) merge(..., all=T), list(all_20k.dd.nna.df, all_220k.dd.nna.df, all_420k.dd.nna.df,
                                                        all_620k.dd.nna.df,all_820k.dd.nna.df,all_1m.dd.nna.df, all_1.2m.dd.nna.df,
                                                        all_1.4m.dd.nna.df, all_1.6m.dd.nna.df,all_1.8m.dd.nna.df))

ggplot(all_merged_finally, aes(x=mean, y=covar)) + theme_bw() + stat_smooth(method=loess,fill="lightblue") + geom_point(aes(color=ReadNum),size=2.5) + xlim(0,5) + ylim(-5,40) + geom_line(x=1,color="green") + geom_line(y=15,color="green") + xlab("Genomic Proportion") + ylab("Coefficient of Variation") + theme(axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16), axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18)) + scale_color_manual(values=c("grey90","grey80","grey70","grey60","grey50","grey40","grey30","grey20","grey10","black"),name="Percent Genome Coverage",breaks=c("20000","220000","420000","620000","820000","1020000","1220000","1420000","1620000","1820000"),labels=c("0.056","0.6","1.2","1.7","2.3","2.8","3.4","3.9","4.5","5.1"))
