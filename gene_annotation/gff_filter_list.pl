#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;
use autodie;
#use Data::Dump::Color;

my (%genehash, %mrnahash);
my $usage  = "$0 list gff";
my $infile = shift or die $usage;
my $gff    = shift or die $usage;

open my $in, '<', $infile;
while (<$in>) {
    chomp;
    my ($gene, $ref, $st, $end) = split /\_/;
    my $mrna;
    ($mrna = $gene) =~ s/gene/mRNA/;
    $genehash{$gene} = 1;
    $mrnahash{$mrna} = 1;
}
close $in;

my $ct = 0;
my $fh;
if ($gff =~ /\.gz$/) {
    open $fh, '-|', 'zcat', $gff;
}
else {
    open $fh, '<', $gff;
}

while (<$fh>) {
    chomp;
    if (/^##\w+/) {
	say $_;
    }
    else {
	my @f = split /\t/;
	if (@f == 9) {
	    my ($gene) = ($f[8] =~ /(gene\d+|mRNA\d+)/);
	    if (defined $gene) {
		if ($gene =~ /gene/) {
		    say "###" if $ct and exists $genehash{$gene};
		    say join "\t", @f and $ct++ if exists $genehash{$gene};
		    
		}
		if ($gene =~ /mRNA/) {
		    say join "\t", @f if exists $mrnahash{$gene};
		}
	    }
	    else {
		say "###";
	    }
	}
    }
}
close $fh;
