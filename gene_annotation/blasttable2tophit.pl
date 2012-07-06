#!/usr/bin/perl -w
#-------------------------------------------------------+
# blasttable2tophit.pl - get only unique hits             |
#-------------------------------------------------------+
# DESCRIPTION: Report only the first (top) hit from a   |
# blasttable (created with -m 8 using legacy BLAST).    |
# This gives the best hit in terms of e-value by        |
# default. To get the best hit by bit score, just sort  |
# the blast report first (i.e. sort -nrk 12).           |
#                                                       |
# AUTHOR: S. Evan Staton                                |
# CONTACT: statonse<at>uga.edu                          |
# STARTED: 2.17.11                                      |
# UPDATED: 2.17.11                                      
#                                                       |
# USAGE:            
#                                                       |
#-------------------------------------------------------+
# TODO: add some simple printing of the stats DONE 2.17.11
#
use strict;

my $usage   = "USAGE: blasttable2tophit.pl blasttablein tophitsout\n\n";
my $infile  = shift or die "\nERROR: No infile was found!\n\n",$usage;
my $outfile = shift or die "\nERROR: No outfile was found!\n\n",$usage;

open( my $in , '<', $infile ) or die "\nERROR: Could not open file: $!\n";
open( my $out, '>' , $outfile ) or die "\nERROR: Could not open file: $!\n";

my $allq = 0;
my $uniqq = 0;
my %seen = ();

while (<$in>) { 
  chomp; 
  $allq++;
  my @blfields = split(/\t/, $_);
  my $key = $blfields[0];
  if (! $seen{ $key }++) {
    print $out $_."\n";
  }
}

$uniqq += keys %seen;

print "$allq Total sequences in report.\n";
print "$uniqq Unique sequences have been written to: $outfile\n";

close($in);
close($out);  

exit;
