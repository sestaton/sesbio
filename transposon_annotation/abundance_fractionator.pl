#!/usr/bin/env perl

## This takes the TE community abundance list that is input to Parthy or Ecolpy
## and returns an integer fraction of the abundance values to make simulations tractable.

use 5.010;
use strict;
use warnings;
use autodie qw(open);
use Getopt::Long;

my $usage = "$0 -i file -n 1000\n";
my $infile;
my $frac;

GetOptions(
	   'i|infile=s' => \$infile,
	   'n|number=i' => \$frac,
	   );

die $usage if !$infile;
$frac //= 10000;

open my $in, '<', $infile;

while (<$in>) {
    chomp;
    my $float   = $_ / $frac;
    my $rounded = int($float + $float/abs($float*2));
    say $rounded;
}
close $in;
