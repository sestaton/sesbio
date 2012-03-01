#!/usr/bin/perl -w

use strict;
#use Data::Dumper;
use Getopt::Long;

my $usage = "USAGE: $0 -i inreport -o outreport

This script takes as input a list of contigs with the length of each
contig separated by tabs such as:

contig-14445    445
contig-8589     589
contig-13776    776
(...)

and the count of each length is output such as (minus the headers):

length  count
1009      1
1038      1
1134      1
1152      1
1236      1
1269      1
130       581
131       554
132       470
133       480
(...)

NB: the lengths will not be reported in the same order as the input file.";

my $infile;
my $outfile;

## initialize counts
#my $total_hits = 0;
#my $unique_hits = 0;

GetOptions(# Required arguments
           'i|infile=s'         => \$infile,
           'o|outfile=s'        => \$outfile,
           );

# open the infile or die with a usage statement
if ($infile && $outfile) {
    open(LENREPORT, '<', $infile) || print "\nERROR: Could not open file: $infile\n";

    open (OUTFILE, '>', $outfile);
}
else {
    if (!$infile){
        die "\nERROR: No infile was given at the command line\n\n",$usage,"\n\n"; 
    }
    if (!$outfile){
        die "\nERROR: No outfile was given at the command line\n\n",$usage,"\n\n";
    }
}

#
# comments must be removed or they will be counted
#
my @contigs = map +(split "\t")[1], <LENREPORT>;
close(LENREPORT);

@contigs = sort { $a <=> $b } @contigs;
my @lengths;
foreach my $len (@contigs) {
    chomp $len;
    push(@lengths,$len);
}
my %seen = ();
my @unique_lengths = grep { ! $seen{$_} ++ } @lengths;   

my $unique = @unique_lengths;
my $total = @contigs;

print "\n","There are: ", $total, " total contigs.\n";
print "\n","There are: ", $unique, " unique contig lengths.\n\n";

count_unique ( @contigs );

close(OUTFILE);

exit;

#
# sub
#
sub count_unique {

    my @array = @_;
    my %count;
    map { $count{$_}++ } @array;

    map {print OUTFILE $_."\t".${count{$_}}."\n"} sort keys(%count);

}


