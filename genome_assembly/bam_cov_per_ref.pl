#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;

# Usage: samtools mpileup file.bam | bam_coverage.pl

my %refs;

while (<>) {
    my @f = split /\t/;
    push @{$refs{$f[0]}}, join "||", @f;
}

for my $ref (keys %refs) {
    my ($len, $num, $min, $max) = ref_stats($refs{$ref});
    printf  "Mean coverage for $ref : %.1f\n", $num/$len;
    printf  "Coverage range for $ref: %d - %d\n", $min, $max;
}

sub ref_stats {
    my ($ref) = @_;

    my $num;     # per residue coverage
    my $len;     # sequence length counter 
    my $min = 0; # minimum coverage
    my $max = 0; # maximum coverage

    for my $aln (@$ref) {
	my ($seq, $pos, $base, $cov, $bases, $quals) = split /\|\|/, $aln;
	$num += $cov;
	$min = $cov if $min > $cov;
	$max = $cov if $max < $cov;
	$len++;
    }

    return ($len, $num, $min, $max);
}
