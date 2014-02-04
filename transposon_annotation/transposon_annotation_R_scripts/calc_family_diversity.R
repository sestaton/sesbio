library(vegan)

fam_bp_div <- read.table("all_15_species_te-families_cval_raw_bp_counts_tab_for_vegan_corr_order_div-format.tsv.txt",sep="\t",row.name=1,header=T)
#H <- diversity(index = "shannon", fam_bp_div)
#S <- specnumber(H)
#J <- H/log(S) # evenness

## to avoid integer overflow errors from R, calculate statistics manually
## these commands below produce the exact same results as those above
hteph.pi <- fam_bp_div[10,]/sum(as.numeric(fam_bp_div[10,]))
hteph.h <- -hteph.pi * log(hteph.pi)
hteph.H <- apply(hteph.h, 1, sum, na.rm = TRUE)
#> hteph.H
#Helianthus niveus ssp. tephrodes 
#                        2.611064 

hteph.J <- hteph.H/log(specnumber(fam_bp_div[10,]))
#> hteph.J
#Helianthus niveus ssp. tephrodes 
#                       0.7031137

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

## RAD and plotting
rad <- radfit(fam_bp_div)
plot(rad)
# change panel order in lattice plot
# plot(fams_rad,index.cond=list(c(4,5,6,13,7,3,14,15,1,12,9,10,11,2,8)))

# or select a single site
#rad4 <- radfit(fam_bp_div[4,])
#plot(rad4)

## calc Arrhenius species-area model and plot quantiles
z <- betadiver(fam_bp_div, "z")
plot(quantile(z))
