library(ggplot2)

## NB: These two plots are completely separate but placed together
##     because the code is obviously very similar

## plot family-level phylogenetic signal
fam_mulitphylostats <- read.table("all_15_species_fam_multiphylosignal_stats_annot_sig-only_corr-order.txt",header=T,sep="\t")
rectangles <- data.frame(xmin = c(0.8,18.8,25.8), xmax = c(18.2, 25.2, 29.2),
                         ymin = c(0, 0, 0), ymax=c(3.5, 3.5, 3.5),
                         t=c("LTR-RT","Non-LTR-RT","Class II"))
dp <- ggplot() + geom_point(data=fam_mulitphylostats, aes(Family, K, color=Superfamily), size=3)
dp + scale_x_discrete(limits=c("COPIA2_I_LC","COPIA2_LTR_MT","Copia15_VV","RLC_X","RLC_amov",
                               "RLC_jiliwu","RLC_ogaow","RLG_X","RLG_kefe","RLG_rewu","RLG_ryse","RLG_teda",
                               "RLG_tewuvu","GYPSY16_I_AG","Gypsy123_I_DR","Gypsy1_SM","DM176","ERV1_N6_I_DR",
                               "L1_11_DR","L1_12_DR","L1_58_ACar","LIN4b_SM","CR1_13_CQ","CR1_58_HM","CR1_79_HM",
                               "P4_AG","SMAR15","ATHPOGON1","Helitron3_PPa"),
                      labels=c("COPIA2_I_LC","COPIA2_LTR_MT","Copia15_VV","RLC_X","RLC_amov","RLC_jiliwu",
                               "RLC_ogaow","RLG_X","RLG_kefe","RLG_rewu","RLG_ryse","RLG_teda","RLG_tewuvu","GYPSY16_I_AG",
                               "Gypsy123_I_DR","Gypsy1_SM","DM176","ERV1_N6_I_DR","L1_11_DR","L1_12_DR","L1_58_ACar",
                               "LIN4b_SM","CR1_13_CQ","CR1_58_HM","CR1_79_HM","P4_AG","SMAR15","ATHPOGON1","Helitron3_PPa")) +
    theme_bw() + theme(axis.text.x=element_text(angle=90, hjust=1, size=8, color="black"), axis.text.y=element_text(color = "black")) + 
    scale_color_manual(name="Superfamily", values=c("Gypsy" = "darkgreen","Copia" = "aquamarine4","L1" = "darkgrey",
                                                    "hAT" = "darkkhaki","Mariner/Tc1" = "orange","L2" = "aquamarine",
                                                    "ERV1" = "darkblue","NeSL" = "chartreuse3","Penelope" = "cyan3",
                                                    "R1" = "black","CR1" = "lightblue","Helitron" = "chocolate")) +
  geom_rect(data=rectangles, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill=t, color="black", alpha=0.2, linetype=0) +
  annotate("text", x=10, y=0.2, size = 3,label="LTR Retrotransposons") +
  annotate("text", x=22, y=0.2, size = 3,label="Non-LTR\nRetrotransposons") +
  annotate("text", x=27.5, y=0.2, size = 3,label="Class II")    

## plot superfamily-level phylogenetic signal
Kstats_reorder <- read.table("all_15_species_superfam_multiphylosignal_stats_sig-only.txt",sep="\t",header=T)
rectangles <- data.frame(xmin = c(0.8,3.8,8.8), xmax = c(3.2, 8.2, 10.2),
                         ymin = c(0, 0, 0), ymax=c(3.5, 3.5, 3.5),
                         t=c("LTR-RT","non-LTR-RT","Class II"))
gp <- ggplot() + geom_point(data=Kstats_reorder, aes(Superfamily, K, color=Superfamily), size=3) +
  theme(axis.title.y = element_blank()) +
  geom_rect(data=rectangles, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill=t, color="black", alpha=0.2, linetype=0)
                          
gp + scale_x_discrete(limits=c("Copia","ERV1","Gypsy","L1","L2","NeSL","Penelope","R1","Mariner/Tc1","hAT"),
                      labels=c("Copia","ERV1","Gypsy","L1","L2","NeSL","Penelope","R1","Mariner/Tc1","hAT")) +
  scale_color_manual(name="Superfamily", values=c("Gypsy" = "darkgreen","Copia" = "aquamarine4","L1" = "darkgrey",
                                           "hAT" = "darkkhaki","Mariner/Tc1" = "chocolate","L2" = "aquamarine",
                                           "ERV1" = "darkblue","NeSL" = "chartreuse3","Penelope" = "cyan3","R1" = "black")) +
  theme_bw() +
  annotate("text", x=2, y=0.2, size = 4,label="LTR Retrotransposons") +
  annotate("text", x=6, y=0.2, size = 4,label="Non-LTR Retrotransposons") +
  annotate("text", x=9.5, y=0.2, size = 4,label="Class II")    
