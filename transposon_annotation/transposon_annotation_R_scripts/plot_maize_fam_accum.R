#library(ggplot2)

maize_cov_sims <- read.table("maize_coverage_stats.tsv",header=T,sep="\t")

ggplot(maize_cov_sims, aes(FamCt, ReadNum)) +
     scale_y_discrete(breaks=c(0,250000,500000,750000,1000000),
                      labels=c("0","250000","500000","750000","1000000")) +
     geom_point() +
     stat_smooth(method="loess",color="red") +
     geom_vline(xintercept = 19.5) +
     geom_hline(yintercept = 350000) +
     geom_rect(aes(xmin = 19.5, xmax = Inf, ymin = 350000, ymax = Inf, fill=T), alpha = 0.010) +
     theme_bw() +
     xlab("TE family number") +
     ylab("Read number") +
     theme(legend.position="none", axis.title.x = element_text(size=14), axis.title.y = element_text(size=14))