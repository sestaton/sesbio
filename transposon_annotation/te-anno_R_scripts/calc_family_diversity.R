library(vegan)

fam_bp_div <- read.table("all_15_species_te-families_cval_raw_bp_counts_tab_for_vegan_corr_order_div-format.tsv.txt",sep="\t",row.name=1,header=T)
H <- diversity(index = "shannon", H)
S <- specnumber(H)
J <- H/log(S) # evenness

## calc Renyi diversities for 6 random sites
k <- sample(nrow(fam_bp_div), 6)
R <- renyi(fam_bp_div[k,])
plot(R)

## generate Fisher's log-series to use as a diversity index
alpha <- fisher.alpha(fam_bp_div)

## calc Fisher's log-series for a plot
fish <- fisherfit(fam_bp_div[1,])
confint(fish) # get confidence intervals

## calc Preston lognormal model
prestondistr(fam_bp_div[4,])

## TE accumulation
teac <- specaccum(fam_bp_div)
#teac # inspect richness and sd
plot(teac, ci.type="polygon", ci.col="yellow")

## get quick diversity stats
specnumber(fam_bp_div)

## Beta diversity
ncol(fam_bp_div)/mean(specnumber(fam_bp_div)) - 1