aster_cval_superfam_dat <- read.table("all_10_species_te-superfamilies_cval_raw_bp_counts_tab.tsv",header=T,sep="\t")
tre_brlen <- read.tree("../aster_wgs_clustering_results/new_annotation_reports_6-17/edited_summaries_6-17/pgls_analysis_with_caper/aster_phylo_with_singles_collapsed_brlen.newick")
aster_cval_superfam_compdat <- comparative.data(tre_brlen, aster_cval_superfam_dat,Species, vcv=TRUE)
aster_cval_pgls_gyp_nonlog <- pgls(formula = log(GenomeSize) ~ log(Gypsy), data = aster_cval_superfam_compdat)