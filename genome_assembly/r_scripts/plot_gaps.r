library(ggplot2)
library(dplyr)

gaps <- read.table("bronze_gap_1mbp_bins.tsv",header=T,sep="\t")

ggplot(gaps %>% mutate(bin.count=1) 
	    %>% group_by(Chr) 
	    %>% mutate(bin.pos = cumsum(bin.count) - 0.5*bin.count), 
	    aes(x=Chr, y=bin.count, fill=N_Perc)) + 
	    coord_flip() + 
	    theme_bw() + 
	    geom_bar(stat="identity") + 
	    scale_x_discrete(limits=c("Ha17","Ha16","Ha15","Ha14","Ha13","Ha12",
	    "Ha11","Ha10","Ha9","Ha8","Ha7","Ha6","Ha5","Ha4","Ha3","Ha2","Ha1")) + 
	    scale_fill_gradient("Percent gaps", low = "white", high = "black", 
	    limits=c(0,max(gaps$N_Perc))) + 
	    xlab("Chromosome") + 
	    ylab("1Mbp bins") + 
	    theme(axis.line = element_blank(),
	    panel.grid.major = element_blank(),
	    panel.grid.minor = element_blank(),
	    panel.border = element_blank(),
	    panel.background = element_blank()) 