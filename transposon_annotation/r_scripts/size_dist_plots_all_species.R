rm(list=ls()) #this clears the workspace to make sure no leftover variables are floating around
graphics.off(); #close all graphics windows

readCls=function(file){
     cls=scan(file,what=character(),sep="\n",comment.char=">",quiet=T) 
     cls=strsplit(cls,split="[ \t]")
     cls
}

# Hann
hann.col <- read.table("Ann1238_cls_size_df_colors.tsv",header=T,sep="\t")
cls=readCls("Ann1238_500k_interl_mgblast_merged_06_17_2013_10_44_12.cls") 
NinClusters=length(unlist(cls))
NinAll=500000
NinSingles=NinAll-NinClusters
clsLength=sapply(cls,length)
clsLength=sort(clsLength,decreasing=T)
clsLength=sapply(cls,length)
clsLength=sort(clsLength,decreasing=T)
clsLengthAll=c(clsLength,rep(1,NinSingles))

png("Ann1238_cls_barplot_color.png")
barplot(clsLength,width=clsLength,space=0,ylim=c(0,max(clsLength)*1.2),col=as.vector(hann.col$Color))
axis(1,at=seq(0,NinAll,length.out=11),label=seq(0,100,by=10));
dev.off()

# Ager
ager.col <- read.table("Ageratina_cls_size_df_colors.tsv",header=T,sep="\t")
cls=readCls("Ageratina_CAGATC_prinseq_trimmed_clean_desc_paired_scr_500k_interl_mgblast_merged_06_17_2013_10_49_12.cls") 
NinClusters=length(unlist(cls))
NinAll=500000
NinSingles=NinAll-NinClusters
clsLength=sapply(cls,length)
clsLength=sort(clsLength,decreasing=T)
clsLength=sapply(cls,length)
clsLength=sort(clsLength,decreasing=T)
clsLengthAll=c(clsLength,rep(1,NinSingles))

png("Ageratina_cls_barplot_color.png")
barplot(clsLength,width=clsLength,space=0,ylim=c(0,max(clsLength)*1.2),col=as.vector(ager.col$Color))
axis(1,at=seq(0,NinAll,length.out=11),label=seq(0,100,by=10));
dev.off()

# Calyc
calyc.col <- read.table("Calyc_size_df_colors.tsv",header=T,sep="\t")
cls=readCls("Calyc_GATCAG_prinseq_trimmed_clean_desc_paired_scr_500k_interl_mgblast_merged_06_17_2013_11_30_16.cls") 
NinClusters=length(unlist(cls))
NinAll=500000
NinSingles=NinAll-NinClusters
clsLength=sapply(cls,length)
clsLength=sort(clsLength,decreasing=T)
clsLength=sapply(cls,length)
clsLength=sort(clsLength,decreasing=T)
clsLengthAll=c(clsLength,rep(1,NinSingles))

png("Calyc_cls_barplot_color.png")
barplot(clsLength,width=clsLength,space=0,ylim=c(0,max(clsLength)*1.2),col=as.vector(calyc.col$Color))
axis(1,at=seq(0,NinAll,length.out=11),label=seq(0,100,by=10));
dev.off()

# CP
cp.col <- read.table("CP_size_df_colors.tsv",header=T,sep="\t")
cls=readCls("CP_AGTCAA_prinseq_trimmed_clean_desc_paired_scr_500k_interl_mgblast_merged_06_17_2013_10_40_34.cls") 
NinClusters=length(unlist(cls))
NinAll=500000
NinSingles=NinAll-NinClusters
clsLength=sapply(cls,length)
clsLength=sort(clsLength,decreasing=T)
clsLength=sapply(cls,length)
clsLength=sort(clsLength,decreasing=T)
clsLengthAll=c(clsLength,rep(1,NinSingles))

png("CP_cls_barplot_color.png")
barplot(clsLength,width=clsLength,space=0,ylim=c(0,max(clsLength)*1.2),col=as.vector(cp.col$Color))
axis(1,at=seq(0,NinAll,length.out=11),label=seq(0,100,by=10));
dev.off()

# Dasy
dasy.col <- read.table("Dasy_size_df_colors.tsv",header=T,sep="\t")
cls=readCls("Dasyphllum_ATCACG_prinseq_trimmed_clean_desc_paired_scr_500k_interl_mgblast_merged_06_17_2013_10_30_53.cls") 
NinClusters=length(unlist(cls))
NinAll=500000
NinSingles=NinAll-NinClusters
clsLength=sapply(cls,length)
clsLength=sort(clsLength,decreasing=T)
clsLength=sapply(cls,length)
clsLength=sort(clsLength,decreasing=T)
clsLengthAll=c(clsLength,rep(1,NinSingles))

png("Dasy_cls_barplot_color.png")
barplot(clsLength,width=clsLength,space=0,ylim=c(0,max(clsLength)*1.2),col=as.vector(dasy.col$Color))
axis(1,at=seq(0,NinAll,length.out=11),label=seq(0,100,by=10));
dev.off()

# Gerb
gerb.col <- read.table("Gerb_size_df_colors.tsv",header=T,sep="\t")
cls=readCls("Gerb_GGCTAC_prinseq_trimmed_clean_desc_paired_scr_500k_interl_mgblast_merged_06_17_2013_10_26_34.cls") 
NinClusters=length(unlist(cls))
NinAll=500000
NinSingles=NinAll-NinClusters
clsLength=sapply(cls,length)
clsLength=sort(clsLength,decreasing=T)
clsLength=sapply(cls,length)
clsLength=sort(clsLength,decreasing=T)
clsLengthAll=c(clsLength,rep(1,NinSingles))

png("Gerb_cls_barplot_color.png")
barplot(clsLength,width=clsLength,space=0,ylim=c(0,max(clsLength)*1.2),col=as.vector(gerb.col$Color))
axis(1,at=seq(0,NinAll,length.out=11),label=seq(0,100,by=10));
dev.off()

# Gnaph
gnaph.col <- read.table("Gnaph_cls_size_df_colors.tsv",header=T,sep="\t")
cls=readCls("Gnaph_ACAGTG_prinseq_trimmed_clean_desc_paired_scr_500k_interl_mgblast_merged_06_17_2013_10_23_13.cls") 
NinClusters=length(unlist(cls))
NinAll=500000
NinSingles=NinAll-NinClusters
clsLength=sapply(cls,length)
clsLength=sort(clsLength,decreasing=T)
clsLength=sapply(cls,length)
clsLength=sort(clsLength,decreasing=T)
clsLengthAll=c(clsLength,rep(1,NinSingles))

png("Gnaph_cls_barplot_color.png")
barplot(clsLength,width=clsLength,space=0,ylim=c(0,max(clsLength)*1.2),col=as.vector(gnaph.col$Color))
axis(1,at=seq(0,NinAll,length.out=11),label=seq(0,100,by=10));
dev.off()

# Saff
saff.col <- read.table("Saff_size_df_colors.tsv",header=T,sep="\t")
cls=readCls("Saff_GATCAG_prinseq_trimmed_clean_desc_paired_scr_500k_interl_mgblast_merged_06_17_2013_10_18_28.cls") 
NinClusters=length(unlist(cls))
NinAll=500000
NinSingles=NinAll-NinClusters
clsLength=sapply(cls,length)
clsLength=sort(clsLength,decreasing=T)
clsLengthAll=c(clsLength,rep(1,NinSingles))

png("Saff_cls_barplot_color.png")
barplot(clsLength,width=clsLength,space=0,ylim=c(0,max(clsLength)*1.2),col=as.vector(saff.col$Color))
axis(1,at=seq(0,NinAll,length.out=11),label=seq(0,100,by=10));
dev.off()

# Sene
sene.col <- read.table("Senecio_size_df_colors.tsv",header=T,sep="\t")
cls=readCls("Senecio_TTAGGC_prinseq_trimmed_clean_desc_paired_scr_500k_interl_mgblast_merged_06_17_2013_09_25_15.cls") 
NinClusters=length(unlist(cls))
NinAll=500000
NinSingles=NinAll-NinClusters
clsLength=sapply(cls,length)
clsLength=sort(clsLength,decreasing=T)
clsLengthAll=c(clsLength,rep(1,NinSingles))

png("Sene_cls_barplot_color.png")
barplot(clsLength,width=clsLength,space=0,ylim=c(0,max(clsLength)*1.2),col=as.vector(sene.col$Color))
axis(1,at=seq(0,NinAll,length.out=11),label=seq(0,100,by=10));
dev.off()

# TKS
tks.col <- read.table("TKS_size_df_colors.tsv",header=T,sep="\t")
cls=readCls("TKS_CAGATC_prinseq_trimmed_clean_desc_paired_scr_500k_interl_mgblast_merged_06_17_2013_09_08_14.cls") 
NinClusters=length(unlist(cls))
NinAll=500000
NinSingles=NinAll-NinClusters
clsLength=sapply(cls,length)
clsLength=sort(clsLength,decreasing=T)
clsLengthAll=c(clsLength,rep(1,NinSingles))

png("TKS_cls_barplot_color.png")
barplot(clsLength,width=clsLength,space=0,ylim=c(0,max(clsLength)*1.2),col=as.vector(tks.col$Color))
axis(1,at=seq(0,NinAll,length.out=11),label=seq(0,100,by=10));
dev.off()
