## If you use this script, please cite: http://www.biomedcentral.com/1471-2164/16/623
#
#
library(ggplot2)

retro_abund <- read.table("retro_abund_genome_size_sort_tab.txt",header=T,sep="\t")

lm_eqn = function(df){
     m = lm(RetroBases ~ GenomeSize, df);
     eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                      list(a = format(coef(m)[1], digits = 2), 
                           b = format(coef(m)[2], digits = 2), 
                           r2 = format(summary(m)$r.squared, digits = 3)))
     as.character(as.expression(eq));                 
}

ggplot(retro_abund, aes(GenomeSize, RetroBases)) +
  geom_point(size=3) + geom_smooth(method="lm",color="red") + 
  theme_bw() +
  xlab("Genome size (Mbp)") +
  ylab("Retrotransposon DNA (Mbp)") +
  annotate("text",y=1200,x=3000,label = lm_eqn(retro_abund), size = 4, parse=TRUE) 