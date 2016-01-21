library(plyr)
> library(ggplot2)
> setwd("Desktop//Hannuus_lines_repeat_analysis")
> lines <- read.table("all_lines_family_stats_6-30.tsv",header=T,sep="\t",comment.char="")
> lines.filtered <- lines[lines$GenomeFrac >= 0.01,]
ggplot(alllines.filt, aes(x=reorder(Identifier, GenomeFrac), y=GenomeFrac)) + geom_bar(aes(fill=Family, order=desc(Family)), stat="identity") + theme_bw() + theme(axis.text.x = element_text(angle=90, hjust=1, color="black"), axis.text.y = element_text(color="black"), axis.title.x = element_blank(), axis.text.y = element_blank())
alllines.noha <- alllines[!alllines$Line == "HA",]
