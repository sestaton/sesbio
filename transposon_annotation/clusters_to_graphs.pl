#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;
use Cwd;
use File::Basename;
use autodie qw(open);

my $usage = "$0 cls_file hitsort cores\n";
my $cls_file = shift or die $usage;
my $hitsort  = shift or die $usage;
my $cores    = shift or die $usage;
#my $gl_dir = shift or die $usage;

clusters2graphs($cls_file, $hitsort, $cores);

exit;
#
# methods
#
sub clusters2graphs {
    my ($cls_file, $hitsort) = @_;

    #unless ($gl_dir =~ /\/$/) {
    #    $gl_dir .= "/";
    #}
    ## open it, write it, then unlink it

    #my $clusters2graph_plot = $cls_file;
    #$clusters2graph_plot .= ".cluster_graph_summary.pdf";
    my $clusters2graph_rscript = $cls_file;
    $clusters2graph_rscript .= ".clusters2graph.rscript";
    my ($hname, $hpath, $hsuffix) = fileparse($hitsort, qr/\.[^.]*/);
    open my $rscript, '>', $clusters2graph_rscript;

    say $rscript "suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(doMC))
suppressPackageStartupMessages(library(multicore))

#########################################################################

## FUNCTIONS :
# this is modification of igraph function  - faster for large graphs
indexing=function(x,ind){
	y=as.numeric(factor(ind,levels=x))
	y
}

graph.data.frame2 <- function(d, directed=TRUE) {
	
	if (ncol(d) < 2) {
		stop(\"the data frame should contain at least two columns\")
	}

	# assign vertex ids
	if (nrow(d)>200000000){
		cat(\'getting unique...\')
		names <- uniqueM(c(d[,1], d[,2]))
		cat(\'done\\n\')

	}else{
		names <- unique( c(as.character(d[,1]), as.character(d[,2])) )

	}
	
	ids <- seq(along=names)

	# create graph
	g <- graph.empty(n=0, directed=directed)
	g <- add.vertices(g, length(ids), name=names)

	# create edge list
	from <- as.character(d[,1])

	to <- as.character(d[,2])
	#edges <- t(matrix(c(ids[from], ids[to]), nc=2))
	edges <- t(matrix(c(indexing(names,from),indexing(names,to)),nc=2))
	
	# edge attributes
	attrs <- list()
	if (ncol(d) > 2) {
		for (i in 3:ncol(d)) {
			newval <- d[,i]
			if (class(newval) == \"factor\") {
				newval <- as.character(newval)
			}
			attrs[[ names(d)[i] ]] <- newval
		}
	}
	
	# add the edges
	g <- add.edges(g, edges, attr=attrs)
	g
}

uniqueM=function(x,n=100)
{ print(length(x)  )
	options(warn=1)
	x=unlist(mclapply(split(x,1:4),FUN=unique))
	options(warn=0)  
	print(length(x))  
	x=unique(x)
	print(length(x)  )
	x
}

#######################
registerDoMC(cores=$cores)

vmax=50000

hitsort=\"$hitsort\"
#clusterlist=\"$cls_file\"
minClsSize=30

cat(\'\\n reading .cls file\\n\')
clusters=scan(file=\"$cls_file\",what=character(),sep=\"\\n\",comment.char=\">\",quiet=T)  # load cluster file
clusters=strsplit(clusters,split=\"[ \\t]\")
clusters=clusters[sapply(clusters,length)>=minClsSize]  # remove smaller clusters
ncolDir=paste(dirname(hitsort),\"/${hname}_ncol\",sep=\'\')

##get hitsort
cat(\"\\nreading hitsort file\\n\\n\")
edgelist=read.table(hitsort,sep=\'\\t\',header=F,as.is=T,colClasses=c(\"character\",\"character\",\"numeric\"))
# create output set directories ncol
dir.create(ncolDir)

#save individual clusters: ncol
for (i in seq_along(clusters)){
     cat(paste(\'saving cluster no.\',i,\'\\n\'))
     subEdgelist=edgelist[(edgelist[,1] %in% clusters[[i]]+edgelist[,2] %in% clusters[[i]])==2,]
     write.table(subEdgelist,sep=\'\\t\',col.names=F,row.names=F,quote=F,file=paste(ncolDir,\"/CL\",i,\'.ncol\',sep=\"\"))
}                                                                                                                                                           
rm(edgelist)           
 
#clusters
# create output set directories GL
GLdir=paste(dirname(hitsort),\"/${hname}_GL\",sep=\'\')
dir.create(GLdir)

# calculate layouts
# read ncol files again:
out=foreach(i = seq_along(clusters),.inorder=FALSE) %dopar% {
	dir.create(ncolDir)
	
	# samples are created in advance!!
	ncolfile=paste(ncolDir,\"/CL\",i,\'.ncol\',sep=\"\")
	ncolSampledFile=paste(ncolDir,\"/CL\",i,\'_sample.ncol\',sep=\"\")
	
	if (file.exists(ncolSampledFile)){
		cat(paste(\'original cluster CL\',i, \'was above threshold!, sample of graph is used\\n\'))
		gd=read.table(file=ncolSampledFile,sep=\'\\t\',header=F,as.is=T,col.names=c(1,2,\'weight\'))
		cat(paste(\'ncol file for cluster CL\',i,\' loaded - \',Sys.time(),\"\\n\",sep=\"\"))
		g=graph.data.frame2(gd,directed=F)
		cat(paste(\'ncol file converted to graph for cluster CL\',i,\" \",Sys.time(),\"\\n\",sep=\"\"))
		rm(gd) # to free memory!
		cat(paste(\'graph of  cluster CL\',i,\' with \',vcount(g),\' nodes, \',ecount(g),\'edges created\',sep=\'\'),\"\\n\")
		CCs=clusters(g)
		#the largest cc
		maxCC=which.max(CCs\$csize)
		g=induced.subgraph(g,vids=which(CCs\$membership==maxCC))
		cat(paste(\'Largest connected component of graph sample of cluster CL\',i,\' - \',vcount(g),\' nodes, \',ecount(g),\'edges created\',sep=\'\'),\"\\n\")
	}else{
		gd=read.table(file=ncolfile,sep=\'\\t\',header=F,as.is=T,col.names=c(1,2,\'weight\'))
		cat(paste(\'ncol file for cluster CL\',i,\' loaded - \',Sys.time(),\"\\n\",sep=\"\"))
		g=graph.data.frame2(gd,directed=F)
		cat(paste(\'ncol file converted to graph for cluster CL\',i,\' \',Sys.time(),\"\\n\",sep=\"\"))
		rm(gd) # to free memory!
		cat(paste(\'graph of  cluster CL\',i,\' with \',vcount(g),\' nodes, \',ecount(g),\'edges created\',sep=\'\'),\"\\n\")
		
	}
			
	cat(paste(\'layout for cluster CL\',i,\' - calculation start at - \',Sys.time(),\"\\n\",sep=\"\"))
	set.seed(1)
	GL=list(G=g,L=layout.fruchterman.reingold(g,dim=3,verbose=F))
	cat(paste(\'layout for cluster CL\',i,\' - calculation finished at - \',Sys.time(),\"\\n\",sep=\"\"))
	save(GL,file=paste(GLdir,\"\/CL\",i,\'.GL\',sep=\"\"))
	cat(paste(\'layout for cluster CL\',i,\' with \',vcount(g),\' nodes, \',ecount(g),\' edges saved at \',Sys.time(),sep=\'\'),\"\\n\")

}";

close $rscript;
	 
system("/usr/local/R/2.15.0/lib64/R/bin/R --vanilla --slave --silent < $clusters2graph_rscript 2> /dev/null");
unlink $clusters2graph_rscript;
}
