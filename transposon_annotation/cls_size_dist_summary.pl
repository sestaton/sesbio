#!/usr/bin/env perl

use strict;
use warnings;
use Cwd;

my $usage = "$0 cls_file seq_file\n";
my $cls_file = shift or die $usage;
my $fas_file = shift or die $usage;

open(my $cls, '<', $cls_file);

my ($seqhash, $seqct) = fas2hash($fas_file);

my $cwd = cwd();

write_cls_size_dist_summary($cls_file, $seqct, $cwd);

#
# subs
#
sub fas2hash {
    my $fas_file = shift;
    open(my $fas, '<', $fas_file);
   
    my $setct = 0;
    local $/ = '>';

    my %seqhash;
    while (my $line = <$fas>) {
        my ($seqid, $seq) = split /\n/, $line;;
        $seqhash{$seqid} = $seq;
	$seqct++ if defined $seq;
    }
    close($fas);
    
    return(\%seqhash, $seqct);
}

sub write_cls_size_dist_summary {
    my ($cls_file, $seqct, $cwd) = @_;

    ## open it, write it, then unlink it

    my $cls_size_dist_plot = $cls_file;
    $cls_size_dist_plot .= ".png";
    my $cls_size_dist_rscript = $cls_file;
    $cls_size_dist_rscript .= ".rscript";
    open(my $rscript, '>', $cls_size_dist_rscript);

    ## 
    ##

    say $rscript "cls=scan(file=\"$cls_file\",what=character(),sep=\"\\n\",comment.char=\">\",quiet=T)
                  cls=strsplit(cls,split=\"[ \\t]\")

png(filename=\"$cls_size_dist_plot\",width=800,height=500)

options(scipen=999) # this turns of printing in scientific notation
#n=${seqct}
#formatC(n, format = \"d\")
#format(n, scientific = FALSE)
#setwd(${cwd})     #################### if this is necessary set it to the cwd
# plot barplot:
NinClusters=length(unlist(cls))
#NinAll=length(seqs)
NinAll=$seqct
NinSingles=NinAll-NinClusters
clsLength=sapply(cls,length)
clsLength=sort(clsLength,decreasing=T)
clsLengthAll=c(clsLength,rep(1,NinSingles))
barplot(clsLengthAll,width=clsLengthAll,space=0,ylim=c(0,max(clsLength)*1.2),ylab=\"number of reads [reads]\",xlab=\"number of reads [%]\",main=paste(NinAll, \"reads total\"))
rect(0,0,NinClusters,clsLength[[1]]*1.2,col=\"#FF000010\")
rect(NinClusters,0,NinAll,clsLength[[1]]*1.2,col=\"#00FF0010\")
text(NinClusters/2,clsLength[[1]]*1.05, labels=paste(NinClusters,\"reads in\\n\",length(cls),\"clusters\"))
text(NinClusters+NinSingles/2,clsLength[[1]]*1.05, labels=paste(NinSingles,\"singlets\"))
axis(1,at=seq(0,NinAll,length.out=11),label=seq(0,100,by=10))

tmp=capture.output(dev.off())";
    close($rscript);

system("R --vanilla --slave --silent < $cls_size_dist_rscript 2> /dev/null");
 
}
