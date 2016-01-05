library(ggplot2)

ha <- read.table("all_lines_family_stats_6-30_ha_over2pc.tsv",comment.char="",header=F,sep="\t")
names(ha) <- c("ID","Line","Oil","Superfamily","Family","PA")

## superfamily summary
ggplot(ha1, aes(x=reorder(ID, PA, sum), y=PA, fill=Superfamily)) + 
	    geom_bar(stat="summary", fun.y=sum,color="black") + 
	    theme_bw() + 
	    theme(axis.text.x=element_blank(), 
	    axis.ticks.x=element_blank(), 
	    axis.title.x=element_text(size=16),
	    axis.title.y=element_text(size=16),
	    axis.text.y=element_text(size=14)) + 
	    ylab("Genome abundance [%]") + 
	    xlab("Line")

## family summary
ggplot(ha1, aes(x=reorder(ID, PA, sum), y=PA, fill=Family)) + 
	    geom_bar(stat="identity",color="black") + 
	    theme_bw() + 
	    theme(axis.text.x=element_blank(), 
	    axis.ticks.x=element_blank(), 
	    axis.title.x=element_text(size=16),
	    axis.title.y=element_text(size=16),
	    axis.text.y=element_text(size=14)) + 
	    ylab("Genome abundance [%]") + 
	    xlab("Line")

