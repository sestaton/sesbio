library(ggplot2)
library(gridExtra)

maize_cov_sims <- read.table("maize_coverage_stats.tsv",header=T,sep="\t")

p <- ggplot(maize_cov_sims, aes(ReadNum, FamCt)) +
     geom_point() +
     stat_smooth(method="loess",color="red") +
     theme_bw() +
     xlab("Read number") +
     ylab("TE family number")

p2 <- ggplot(maize_cov_sims, aes(ReadNum, RepeatTotFromAnnot)) +
    scale_x_continuous(breaks=c(0,250000,500000,750000,1000000),labels=c("0","250,000 (1.0)","500,000 (2.0)","750,000 (3.0)","1,000,000 (4.0)")) +
    geom_point() +
    stat_smooth(method="loess",color="red") +
    theme_bw() +
    xlab("Read number (percent genome coverage)") +
    ylab("TE fraction of genome")

gA <- ggplot_gtable(ggplot_build(p + theme(legend.position="none",
                                           axis.text.x = element_blank(),
					   axis.title.x = element_blank(),
					   axis.ticks.x = element_blank(),
					   axis.text.y = element_text(size=12),
					   axis.title.y = element_text(size=14))))
gB <- ggplot_gtable(ggplot_build(p2 + theme(legend.position="none",
					    axis.text.x = element_text(size=12),
					    axis.text.y = element_text(size=12),
					    axis.title.y = element_text(size=14),
                                            axis.title.x = element_text(color = "black", size = 14))))

maxWidth = grid::unit.pmax(gA$widths[2:3], gB$widths[2:3])
gA$widths[2:3] <- as.list(maxWidth)
gB$widths[2:3] <- as.list(maxWidth)
print(grid.arrange(gA, gB, ncol=1))
