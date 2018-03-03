#!/usr/bin/env perl

# blasttable2tophit.pl - get only unique hits

# DESCRIPTION: Report only the first (top) hit from a
# blasttable (created with -m 8 using legacy BLAST).
# This gives the best hit in terms of e-value by    
# default. To get the best hit by bit score, just sort
# the blast report first (i.e. sort -nrk 12).  
#
use 5.010;
use strict;
use warnings;

my $usage   = "USAGE: blasttable2tophit.pl blasttablein tophitsout\n\n";
my $infile  = shift or die "\nERROR: No infile was found!\n\n",$usage;
my $outfile = shift or die "\nERROR: No outfile was found!\n\n",$usage;

open my $in , '<', $infile or die "\nERROR: Could not open file: $!\n";
open my $out, '>' , $outfile or die "\nERROR: Could not open file: $!\n";

my $allq = 0;
my $uniqq = 0;
my %seen = ();

while (my $line = <$in>) { 
  chomp $line; 
  $allq++;
  my @blfields = split /\t/, $line;
  my $key = $blfields[0];
  if (! $seen{ $key }++) {
    say $out $line;
  }
}

$uniqq += keys %seen;

say "$allq Total sequences in report.";
say "$uniqq Unique sequences have been written to: $outfile";

close $in;
close $out;  

exit;
