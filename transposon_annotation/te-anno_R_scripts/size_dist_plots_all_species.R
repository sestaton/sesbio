rm(list=ls()) #this clears the workspace to make sure no leftover variables are floating around
graphics.off(); #close all graphics windows

readCls=function(file){
     cls=scan(file,what=character(),sep="\n",comment.char=">",quiet=T) 
     cls=strsplit(cls,split="[ \t]")
     cls
}

# Hann
hann.col <- read.table("Ann1238_1m_interl_mgblast_anno_perc_df_colors.tsv",header=T,sep="\t")
cls=readCls("Ann1238_1m_interl_mgblast_merged_04_21_2013_10_37_46.cls") 
NinClusters=length(unlist(cls))
NinAll=1000000
NinSingles=NinAll-NinClusters
clsLength=sapply(cls,length)
clsLength=sort(clsLength,decreasing=T)

png("Ann1238_cls_barplot_color.png")
clsLength=sapply(cls,length)
clsLength=sort(clsLength,decreasing=T)
clsLengthAll=c(clsLength,rep(1,NinSingles))
barplot(clsLength,width=clsLength,space=0,ylim=c(0,max(clsLength)*1.2),col=as.vector(hann.col$Color))
axis(1,at=seq(0,NinAll,length.out=11),label=seq(0,100,by=10));

# Ager
ager.col <- read.table("Ageratina_anno_perc_df_colors.tsv",header=T,sep="\t")
cls=readCls("Ageratina_CAGATC_prinseq_trimmed_clean_desc_paired_scr_1m_interl_mgblast_merged_04_21_2013_12_24_23.cls") 
NinClusters=length(unlist(cls))
NinAll=1000000
NinSingles=NinAll-NinClusters
clsLength=sapply(cls,length)
clsLength=sort(clsLength,decreasing=T)

png("Ageratina_cls_barplot_color.png")
clsLength=sapply(cls,length)
clsLength=sort(clsLength,decreasing=T)
clsLengthAll=c(clsLength,rep(1,NinSingles))
barplot(clsLength,width=clsLength,space=0,ylim=c(0,max(clsLength)*1.2),col=as.vector(ager.col$Color))
axis(1,at=seq(0,NinAll,length.out=11),label=seq(0,100,by=10));

# Calyc
calyc.col <- read.table("Calyc_anno_perc_df_colors.tsv",header=T,sep="\t")
cls=readCls("Calyc_500k_interl_mgblast_merged_04_22_2013_17_34_22.cls") 
NinClusters=length(unlist(cls))
NinAll=500000
NinSingles=NinAll-NinClusters
clsLength=sapply(cls,length)
clsLength=sort(clsLength,decreasing=T)

png("Calyc_anno_perc_df_colors.tsv")
clsLength=sapply(cls,length)
clsLength=sort(clsLength,decreasing=T)
clsLengthAll=c(clsLength,rep(1,NinSingles))
barplot(clsLength,width=clsLength,space=0,ylim=c(0,max(clsLength)*1.2),col=as.vector(calyc.col$Color))
axis(1,at=seq(0,NinAll,length.out=11),label=seq(0,100,by=10));

# CP
cp.col <- read.table("CP_anno_perc_df_colors.tsv",header=T,sep="\t")
cls=readCls("CP_AGTCAA_prinseq_trimmed_clean_desc_paired_scr_1m_interl_mgblast_merged_04_20_2013_13_59_55.cls") 
NinClusters=length(unlist(cls))
NinAll=1000000
NinSingles=NinAll-NinClusters
clsLength=sapply(cls,length)
clsLength=sort(clsLength,decreasing=T)

png("CP_cls_barplot_color.png")
clsLength=sapply(cls,length)
clsLength=sort(clsLength,decreasing=T)
clsLengthAll=c(clsLength,rep(1,NinSingles))
barplot(clsLength,width=clsLength,space=0,ylim=c(0,max(clsLength)*1.2),col=as.vector(cp.col$Color))
axis(1,at=seq(0,NinAll,length.out=11),label=seq(0,100,by=10));

# Dasy
dasy.col <- read.table("Dasy_anno_perc_df_colors.tsv",header=T,sep="\t")
cls=readCls("Dasyphllum_ATCACG_prinseq_trimmed_clean_desc_paired_scr_1m_interl_mgblast_merged_04_21_2013_02_30_50.cls") 
NinClusters=length(unlist(cls))
NinAll=1000000
NinSingles=NinAll-NinClusters
clsLength=sapply(cls,length)
clsLength=sort(clsLength,decreasing=T)

png("Dasy_cls_barplot_color.png")
clsLength=sapply(cls,length)
clsLength=sort(clsLength,decreasing=T)
clsLengthAll=c(clsLength,rep(1,NinSingles))
barplot(clsLength,width=clsLength,space=0,ylim=c(0,max(clsLength)*1.2),col=as.vector(dasy.col$Color))
axis(1,at=seq(0,NinAll,length.out=11),label=seq(0,100,by=10));

# Gerb
gerb.col <- read.table("Gerb_anno_perc_df_colors.tsv",header=T,sep="\t")
cls=readCls("Gerb_GGCTAC_prinseq_trimmed_clean_desc_paired_scr_1m_interl_mgblast_merged_04_21_2013_04_44_23.cls") 
NinClusters=length(unlist(cls))
NinAll=1000000
NinSingles=NinAll-NinClusters
clsLength=sapply(cls,length)
clsLength=sort(clsLength,decreasing=T)

png("Gerb_cls_barplot_color.png")
clsLength=sapply(cls,length)
clsLength=sort(clsLength,decreasing=T)
clsLengthAll=c(clsLength,rep(1,NinSingles))
barplot(clsLength,width=clsLength,space=0,ylim=c(0,max(clsLength)*1.2),col=as.vector(gerb.col$Color))
axis(1,at=seq(0,NinAll,length.out=11),label=seq(0,100,by=10));

# Gnaph
gnaph.col <- read.table("Gnaph_anno_perc_df_colors.tsv",header=T,sep="\t")
cls=readCls("Gnaph_ACAGTG_prinseq_trimmed_clean_desc_paired_scr_1m_interl_mgblast_merged_04_21_2013_11_27_12.cls") 
NinClusters=length(unlist(cls))
NinAll=1000000
NinSingles=NinAll-NinClusters
clsLength=sapply(cls,length)
clsLength=sort(clsLength,decreasing=T)

png("Gnaph_cls_barplot_color.png")
clsLength=sapply(cls,length)
clsLength=sort(clsLength,decreasing=T)
clsLengthAll=c(clsLength,rep(1,NinSingles))
barplot(clsLength,width=clsLength,space=0,ylim=c(0,max(clsLength)*1.2),col=as.vector(gnaph.col$Color))
axis(1,at=seq(0,NinAll,length.out=11),label=seq(0,100,by=10));

# Saff
saff.col <- read.table("Saff_anno_perc_df_colors.tsv",header=T,sep="\t")
cls=readCls("Saff_GATCAG_prinseq_trimmed_clean_desc_paired_scr_1m_interl_mgblast_merged_04_20_2013_14_00_21.cls") 
NinClusters=length(unlist(cls))
NinAll=1000000
NinSingles=NinAll-NinClusters
clsLength=sapply(cls,length)
clsLength=sort(clsLength,decreasing=T)

png("Saff_cls_barplot_color.png")
clsLength=sapply(cls,length)
clsLength=sort(clsLength,decreasing=T)
clsLengthAll=c(clsLength,rep(1,NinSingles))
barplot(clsLength,width=clsLength,space=0,ylim=c(0,max(clsLength)*1.2),col=as.vector(saff.col$Color))
axis(1,at=seq(0,NinAll,length.out=11),label=seq(0,100,by=10));

# Sene
sene.col <- read.table("Sene_anno_perc_df_colors.tsv",header=T,sep="\t")
cls=readCls("Senecio_TTAGGC_prinseq_trimmed_clean_desc_paired_scr_1m_interl_mgblast_merged_04_20_2013_13_59_41.cls") 
NinClusters=length(unlist(cls))
NinAll=1000000
NinSingles=NinAll-NinClusters
clsLength=sapply(cls,length)
clsLength=sort(clsLength,decreasing=T)

png("Sene_cls_barplot_color.png")
clsLength=sapply(cls,length)
clsLength=sort(clsLength,decreasing=T)
clsLengthAll=c(clsLength,rep(1,NinSingles))
barplot(clsLength,width=clsLength,space=0,ylim=c(0,max(clsLength)*1.2),col=as.vector(sene.col$Color))
axis(1,at=seq(0,NinAll,length.out=11),label=seq(0,100,by=10));

# TKS
tks.col <- read.table("TKS_anno_perc_df_colors.tsv",header=T,sep="\t")
cls=readCls("TKS_CAGATC_prinseq_trimmed_clean_desc_paired_scr_1m_interl_mgblast_merged_04_20_2013_19_31_05.cls") 
NinClusters=length(unlist(cls))
NinAll=1000000
NinSingles=NinAll-NinClusters
clsLength=sapply(cls,length)
clsLength=sort(clsLength,decreasing=T)

png("TKS_cls_barplot_color.png")
clsLength=sapply(cls,length)
clsLength=sort(clsLength,decreasing=T)
clsLengthAll=c(clsLength,rep(1,NinSingles))
barplot(clsLength,width=clsLength,space=0,ylim=c(0,max(clsLength)*1.2),col=as.vector(tks.col$Color))
axis(1,at=seq(0,NinAll,length.out=11),label=seq(0,100,by=10));

