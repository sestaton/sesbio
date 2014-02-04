library(ggplot2)
library(gridExtra)

# modified from: http://stackoverflow.com/a/7549819/1543853
lm_eqn = function(df){
     m = lm(Families ~ Genome_size, df);
     eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                      list(a = format(coef(m)[1], digits = 2), 
                           b = format(coef(m)[2], digits = 2), 
                           r2 = format(summary(m)$r.squared, digits = 3)))
     as.character(as.expression(eq));                 
}

lm2_eqn = function(df){
     m = lm(Family_mean ~ Genome_size, df);
     eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,
                      list(a = format(coef(m)[1], digits = 2),
                           b = format(coef(m)[2], digits = 2),
                           r2 = format(summary(m)$r.squared, digits = 3)))
     as.character(as.expression(eq));
}

fam_size_cval_reps <- read.table("Aster_15_species_famct_cval_fammean_all_reps.txt",header=T,sep="\t")

fit2 <- lm(Families ~ Genome_size, data = fam_size_cval_reps)

p <- ggplot(fit2$model, aes_string(x = names(fit2$model)[2], y = names(fit2$model)[1])) + geom_point(size=3) + stat_smooth(method="lm", col="red") + 
  annotate("text",x= 3500, y=38, label = lm_eqn(fam_size_cval_reps), size = 3, parse=TRUE)

fit3 <- lm(Family_mean ~ Genome_size, data = fam_size_cval_reps)
p2 <- ggplot(fit3$model, aes_string(x = names(fit3$model)[2], y = names(fit3$model)[1])) + geom_point(size=3) + stat_smooth(method="lm", col="red") +
   annotate("text",x= 3500, y=0.95, label = lm2_eqn(fam_size_cval_reps), size = 3, parse=TRUE)

## align plots
gA <- ggplot_gtable(ggplot_build(p + theme_bw() + ylab("TE family richness") + theme(axis.title.x = element_blank())))
gB <- ggplot_gtable(ggplot_build(p2 + theme_bw() + ylab("Mean TE family size (Genome %)") + xlab("Genome size (Mbp)")))

maxWidth = grid::unit.pmax(gA$widths[2:3], gB$widths[2:3])
gA$widths[2:3] <- as.list(maxWidth)
gB$widths[2:3] <- as.list(maxWidth)
print(grid.arrange(gA, gB, ncol=1))