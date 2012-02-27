#!/usr/bin/perl -w

use strict;
#use Data::Dumper;
use Getopt::Long;

my $usage = "USAGE: get_unique_word_counts.pl -i inreport -o outreport

This script takes as input a list of words that may be separated by spaces;
the words are evaluated by unique lines.";

my $infile;
my $outfile;

## initialize counts
#my $total_hits = 0;
#my $unique_hits = 0;

GetOptions(# Required arguments
           "i|infile=s"         => \$infile,
           "o|outfile=s"        => \$outfile,
           );

# open the infile or die with a usage statement
if ($infile && $outfile) {
    open(BLASTREPORT, $infile)|| print "Error: can't read $infile\n";

    open (OUTFILE, ">$outfile");
}
else {
    if (!$infile){
        die "\n","ERROR: No infile was given at the command line\n\n",$usage,"\n\n"; 
    }
    if (!$outfile){
        die "\n","ERROR: No outfile was given at the command line\n\n",$usage,"\n\n";
    }
}

#
# comments must be removed or they will be counted
#
my @repnames = map +(split "\n")[0], <BLASTREPORT>;

my %seen = ();
my @unique_repnames = grep { ! $seen{$_} ++ } @repnames;   # preserves the order of elements
close(BLASTREPORT);

my $unique = @unique_repnames;
my $query = @repnames;

print "\n","There are: ", $query, " total repbase blast hits\n";
print "\n","There are: ", $unique, " unique repbase blast hits\n\n";

count_unique ( @repnames );

sub count_unique {

    my @array = @_;
    my %count;
    map { $count{$_}++ } @array;

      #print

    map {print OUTFILE $_."\t".${count{$_}}."\n"} sort keys(%count);

}

close(OUTFILE);

exit;
