library("ggplot2")

sims <- read.table("all_te_similarity_0106.txt",sep="\t",header=F)
names(sims) <- c("type","element","length","similarity")
sims$type <- factor(sims$type, 
	     	  levels = c("copia","gypsy","unclassified-ltr","trim","hAT",
		  "mutator","tc1-mariner","unclassified-tir"), 
		  labels = c("Copia","Gypsy","Unclassified-LTR","TRIM","hAT",
		  "Mutator","Tc1-Mariner","Unclassified-TIR"))

ggplot(sims, aes(similarity, y=..density..)) + 
	     geom_histogram(binwidth=.3,fill="cornsilk",colour="grey60",size=.2) + 
	     geom_density(size=1,alpha=.3) + 
	     scale_x_reverse() + 
	     theme_bw() + 
	     facet_grid(type ~ .) + 
	     ylab("Density") + 
	     xlab("LTR/TIR similarity")