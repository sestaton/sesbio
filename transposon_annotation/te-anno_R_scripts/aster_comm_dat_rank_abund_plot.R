library(vegan)
superfam_bp <- read.table("all_superfams_bp_tab_ex.txt",sep="\t",row.names=1,header=T)
superfam_bp_rad <- radfit(superfam_bp)
plot(superfam_bp_rad)