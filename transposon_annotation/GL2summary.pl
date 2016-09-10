#!/usr/bin/env perl

## R code below is from the SeqGrapheR package by Petr Novak
## http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-378

use 5.010;
use strict;
use warnings;
use autodie qw(open);
use Cwd;

my $usage    = "$0 cls_file gl_dir\n";
my $cls_file = shift or die $usage;
my $gl_dir   = shift or die $usage;

GL2summary($gl_dir, $cls_file);

# methods
sub GL2summary {
    my ($gl_dir, $cls_file) = @_;

    unless ($gl_dir =~ /\/$/) {
        $gl_dir .= "/";
    }

    ## edit these filenames to something better
    my $gl2summary_plot = $cls_file;
    $gl2summary_plot .= ".cluster_graph_summary.pdf";
    my $gl2summary_rscript = $cls_file;
    $gl2summary_rscript .= "gl2summary.rscript";
    open my $rscript, '>', $gl2summary_rscript;

    say "gl2summary_rscript $gl2summary_rscript say glsummary_plot $gl2summary_plot";
    
say $rscript "suppressPackageStartupMessages(library(igraph))
plotg=function(GG,LL,wlim=NULL,...){
	
	e=get.edgelist(GG,names=F)
	w=E(GG)\$weight
	if (!is.null(wlim)) {e=e[w>wlim,]; w=w[w>wlim]}
	X0=LL[e[,1],1]
	Y0=LL[e[,1],2]
	X1=LL[e[,2],1]
	Y1=LL[e[,2],2]
	plot(range(LL[,1]),range(LL[,2]),xlab=\"\",ylab=\"\",axes=F,type=\"n\",...)
	brv=\'grey\'
	segments(X0,Y0,X1,Y1,lwd=.5,col=brv)
	points(LL,pch=18,cex=.4,...)
}
############################################################################
Diameter=''
tmp=capture.output(as.null(Diameter)) # we don't want to calculate this because it takes too long
# get .GL file names
setwd(\"$gl_dir\")
GLfiles=system(\"ls *.GL\",intern=T)
# assume it is numbered CLXX.GL - get number:
GLfilesIndex=as.numeric(gsub(\"[^0-9]\",\"\",GLfiles))
GLfiles=GLfiles[order(GLfilesIndex)]
			
j=0
page=1
figPerPage=12
ll=c(4,4)
pngFiles=list()
newPage=seq(1,length(GLfiles),figPerPage)
pageLayout=matrix(c(matrix(1:8,ncol=2,byrow=T),matrix(9:16,ncol=2,byrow=T),matrix(17:24,ncol=2,byrow=T)),ncol=4,byrow=T)
	
for (i in GLfiles){
	j=j+1
	if (j %in% newPage){
		graphics.off()
		pngFiles[[j]]=paste(\'page\',sprintf(\"%04d\", page),\".png\",sep=\'\')
		png(pngFiles[[j]],width=2481,height=3507,pointsize=40)
		layout(pageLayout,heights=c(2,1,2,1,2,1,2,1))
		par(mar=c(0.5,0.5,0.5,0.5))
		page=page+1
	}
	cat (i,\'...\')
	load(i)
	plotg(GL\$G,GL\$L)
	par(mar=c(0.0,0.0,0.0,0.0))
	plot(0:10,0:10,type=\'n\',axes=F,xlab=\'\',ylab=\'\',main=paste(\'\\nCL\',j,sep=\'\'))
	text2plot=paste(\'Number of reads:   \',vcount(GL\$G),\'\\nNumber of pairs:   \',ecount(GL\$G),
			\"\\nDensity:      \",signif(graph.density(GL\$G),4),\"\\nDiameter:      \",ifelse(is.null(Diameter),\"NA\",length(get.diameter(GL\$G,directed=F))),
       			\"\\nMean edge weigth:   \", signif(mean(E(GL\$G)\$weight),5),\"\\nMax. degree:     \",max(degree(GL\$G)),sep=\'\')
       	text(1,4,labels=text2plot,pos=4)
       	abline(h=0)
       	par(mar=c(0.5,0.5,0.5,0.5))
       	cat(\"done\\n\")
}
tmp2=capture.output(dev.off())
# make pdf file from pngs 
warnings()
cmd=paste(\"convert page????.png $gl2summary_plot\")
system(cmd)
tmp=lapply(pngFiles,unlink)  # remove all png files";

    close $rscript;

    system("R --vanilla --slave --silent < $gl2summary_rscript 2> /dev/null");
    #unlink($rscript); ## keep for debugging
}
