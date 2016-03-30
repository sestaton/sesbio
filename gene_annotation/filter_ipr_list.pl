#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;
use autodie;
#use Data::Dump::Color;

my (%genehash, %mrnahash);
my $usage  = "$0 list ipr";
my $infile = shift or die $usage;
my $ipr    = shift or die $usage;

open my $in, '<', $infile;
while (<$in>) {
    chomp;
    my ($gene, $ref, $st, $end) = split /\_/;
    my $mrna;
    $genehash{$gene} = 1;
}
close $in;

open my $fh, '<', $ipr;
while (<$fh>) {
    chomp;
    my @f = split /\t/;
    if (exists $genehash{$f[0]}) {
	say join "\t", @f;
    }
}
close $fh;
