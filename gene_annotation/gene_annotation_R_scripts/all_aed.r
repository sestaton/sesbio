library(ggplot2)
aed <- read.table("all_aed_scores.txt",header=F,sep=" ")
names(aed) <- c("sample","AED")

ggplot(aed, aes(x = AED, color=sample)) +
  stat_ecdf(aes(x = AED), geom="line",lwd=1.4) +
  xlim(0,1.0) +
  theme_bw() +
  xlab("AED") +
  ylab("Cumulative fraction of gene annotations") +
  scale_color_discrete(limits=c("makerre","makerst","sunflower","tair10"),
                   labels=c("MAKER Arab. update","MAKER Arab. de novo",
                            "MAKER sunflower","TAIR10 annotations")) +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 12),
        legend.title = element_blank())