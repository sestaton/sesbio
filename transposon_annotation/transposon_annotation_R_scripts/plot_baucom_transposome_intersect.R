library(ggplot2)

zeamelt <- read.table("maize_transposome_1m_5-4_top20_melt_tab.txt",header=T)

gp <- ggplot(zeamelt, aes(Family, GenomePerc*100, color=Study, group=1)) + geom_point(size=3) + geom_smooth(color="red") + theme_bw() 
gp + ylab("Genome percent [%]") +
  scale_color_manual(name="Study", values = c("Baucom" = "black", "Staton" = "blue")) +
   scale_x_discrete(limits=c("RLG_huck","RLC_ji","RLG_cinful","RLC_opie","RLG_flip","RLG_xilon","RLG_prem",
                      "RLG_gyma","RLG_grande","RLG_doke","RLC_giepum","RLG_puck",
                      "RLG_tekay","RLG_uwum","RLG_dagaf","RLC_wiwa","RLG_CRM"),
                    labels=c("RLG_huck","RLC_ji","RLG_cinful","RLC_opie","RLG_flip","RLG_xilon","RLG_prem",
                      "RLG_gyma","RLG_grande","RLG_doke","RLC_giepum","RLG_puck",
                      "RLG_tekay","RLG_uwum","RLG_dagaf","RLC_wiwa","RLG_CRM")) +
  theme(axis.text.x=element_text(angle=90, hjust=1, size=14, color="black"),
        axis.text.y=element_text(size=14, color = "black"),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=16, color = "black"))
