#!/usr/bin/env perl

## This takes the TE community abundance list that is input to Parthy or Ecolpy
## and returns an integer fraction of the abundance values to make simulations tractable.

use 5.010;
use strict;
use warnings;

my $frac = shift;

$frac //= 10000;

while (<>) {
    chomp;
    my $float   = $_ / $frac;
    my $rounded = int($float + $float/abs($float*2));
    say $rounded;
}
