helianthus_labeller <- function(variable,value) {
  return(helianthus_names[value])
}

hel_fam_abund <- read.table("all_helianthus_and_phoeb_families_sorted_by_abund_reduced_names.txt",header=T,sep="\t")

## order panels
dat <- within(hel_fam_abund,
              Species <- factor(Species, levels =c("Helianthus_argophyllus","Helianthus_annuus",
                                           "Helianthus_niveus_ssp._tephrodes","Helianthus_verticillatus",
                                           "Helianthus_porteri","Phoebanthus_tenuifolius")))

## inspect levels
# > unique(levels(dat$Family))

# order labels in legend
dat$Family <- factor(dat$Family, levels=c("RLC_amov","RLC_jiliwu","RLC_ogaow","RLG_begi",
                                   "RLG_iketas","RLG_kefe","RLG_rahi","RLG_rewu",
                                   "RLG_teda","RLG_wily","RLG_wimu"))

# create custom labels for facet titles
helianthus_names <- list(
                         "Helianthus_argophyllus"="H. argophyllus",
                         "Helianthus_annuus"="H. annuus",
                         "Helianthus_niveus_ssp._tephrodes"="H. niveus ssp. tephrodes",
                         "Helianthus_verticillatus"="H. verticillatus",
                         "Helianthus_porteri"="H. porteri",
                         "Phoebanthus_tenuifolius"="Phoebanthus tenuifolius"
                         )

# the reorder sorts by total %, not within each facet
ggplot(dat, aes(x=GPerc, y=reorder(Family, GPerc))) +
    theme_bw() +
    geom_segment(aes(yend=Family), xend=0, colour="grey50") +
    geom_point(pch=21, size=3, aes(fill=Family)) +
    scale_fill_brewer(palette="Paired") +
    facet_grid(Species ~ ., labeller=helianthus_labeller, scales="free_y", space="free_y") +
    theme(strip.text.x = element_text(size = 6, colour = "black")) +
    xlab("Mean genomic abundance [%]") +
    ylab("Family")
