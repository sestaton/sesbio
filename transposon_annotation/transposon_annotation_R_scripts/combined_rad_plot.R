## Fitting each of these models was only possible by changing the
## code for each model in vegan to specify we are summing integers.
## Otherwise, we get integer overflow errors when trying to sum
## the data.

library(vegan)

fam_bp_div_comb <- read.table("aster_15_species_combined_comparative_df_for_caper_families_tab.txt",header=T,sep="\t",row.name=1)

## fit all models
fam_bp_div_comb_rad <- radfit(fam_bp_div_comb)

## save the data
save(fam_bp_div_comb_rad, file = "modified_vegan_radfit_all15_species.RData")

## load the data (in RStudio for fine-tuning plots)
#load("modified_vegan_radfit_all15_species.RData")

## correct order phylogenetically
#plot(fam_bp_div_comb_rad,index.cond=list(c(8,2,10,11,9,12,1,7,14,15,3,13,6,5,4)))