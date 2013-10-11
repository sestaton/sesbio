# CL_repeat_plot.0.4.R
# cluster annotations
rm(list=ls()) #this clears the workspace to make sure no leftover variables are floating around
graphics.off(); #close all graphics windows

CLannot <- read.table("/Users/spencerstaton/Desktop/eR/misc\ scripts/ha412ho-500k-split1_self_mgblast_out_parsed-cluster_size_report.txt",header=F,sep="\t")

cluster <- CLannot[,1]
readcount <- CLannot[,2]
All <- length(readcount)

tot      <- sum(readcount)
mincount <- min(readcount)
maxcount <- max(readcount)

percentg <- readcount/tot * 100
cumperc  <- sum(percentg)
teperc   <- cumperc/4
#cat("here is readcount",readcount)
#cat("Contents of percentg",percentg,fill=TRUE)
#cat("Cumulative perct",cumperc,fill=TRUE)
#cat("Total read count for all clusters: ",tot,fill=TRUE)
#brk <- readcount

#-----------------------------------------
# from 454Summary.R
# plot barplot:
NinClusters=457845                  # number in clusters (457845)
NinAll=540574                       # number in total (540574)
NinSingles=NinAll-NinClusters
#clsLength=sapply(readcount,length)
#clsLength=sort(clsLength,decreasing=T)
clsLengthAll=c(readcount,rep(1,NinSingles))

Totlen=sum(clsLengthAll)
PerTot=clsLengthAll/Totlen * 100

#SingLen=sum(NinSingles)
SingPerTot=NinSingles/Totlen * 100

SumPercTot=sum(PerTot)
MultiCopy=SumPercTot-SingPerTot
#cat("Vector of percents is: ",PerTot)

# Calculate vector of percents
NinPer=NinAll*0.90
round(NinPer)
#formatC(NinPer, mode="double", format="d")
EigPer=NinAll*0.80
SevPer=NinAll*0.70
SixPer=NinAll*0.60
FivPer=NinAll*0.50
FouPer=NinAll*0.40
ThrPer=NinAll*0.30
TwoPer=NinAll*0.20
TenPer=NinAll*0.10
#options(digits=0)
repeatcolors=c("darkgreen","darkgreen","darkgreen","darkgreen","darkgreen","darkgreen","aquamarine4","darkgreen","darkgreen","darkgreen","aquamarine4","aquamarine4","darkgreen","darkgreen","darkgreen","darkolivegreen3","darkgreen","aquamarine4","darkgreen",
               "aquamarine4","chartreuse","darkgreen","aquamarine4","darkgreen","aquamarine4","aquamarine4","aquamarine4","aquamarine4","darkgreen","aquamarine4","aquamarine","darkgreen","azure2","darkgreen","darkgreen","aquamarine4","chartreuse","darkgreen","aquamarine4",
               "darkgreen","darkgreen","aquamarine4","darkgreen","aquamarine4","chartreuse","chartreuse","chartreuse","chartreuse","darkgreen","darkgreen")

par(mar=c(5, 4, 4, 4) + 0.1)

barplot(PerTot,width=PerTot,ylim=c(0,max(PerTot)*1.2),space=0,col=repeatcolors,xaxt="n",ylab="Percentage of all reads [%]",xlab="Cumulative percentage of all reads [%]")

#text(NinClusters/2,clsLength[[1]]*1.05, labels=paste(NinClusters,"reads in\n",length(cls),"clusters")) 
#text(NinClusters+NinSingles/2,clsLength[[1]]*1.05, labels=paste(NinSingles,"singlets"))

#rect(0,0,SumPercTot,max(PerTot)*1.2,col="#FF000010")
#rect(MultiCopy,0,SumPercTot,max(PerTot)*1.2,col="#00FF0010")
rect(0,0,SumPercTot,max(PerTot)*1.2)
rect(MultiCopy,0,SumPercTot,max(PerTot)*1.2)

axis(1,at=seq(0,SumPercTot,length.out=11),label=seq(0,100,by=10))
#axis(3,at=seq(10,SumPercTot-10,by=10),label=c(54057,108114,162172,216229,FivPer,324344,378401,432459,486516))
axis(4,at=seq(0,8,by=2),label=c(0,15000,25000,35000,45000))
#mtext("Cumulative count of reads [n]",side=3,line=2.5)
mtext("Count of reads [n]",side=4,line=2.5)
text(92,max(PerTot)*1.05, labels=paste(NinSingles,"singletons"))
legend(60,8,c("Ty3/gypsy","Ty1/copia","Non-LTR","Helitron","Class II","rRNA"),lwd=4,col=c("darkgreen","aquamarine4","aquamarine","darkolivegreen3","chartreuse","azure2"),bty="n")
#box()
