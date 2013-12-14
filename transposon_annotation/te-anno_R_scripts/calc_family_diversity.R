library(vegan)

fam_bp_div <- read.table("all_15_species_te-families_cval_raw_bp_counts_tab_for_vegan_corr_order_div-format.tsv.txt",sep="\t",row.name=1,header=T)
H <- diversity(index = "shannon", H)
S <- specnumber(H)
J <- H/log(S) # evenness