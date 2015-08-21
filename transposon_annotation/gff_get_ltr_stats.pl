#!/usr/bin/env perl

use strict;
use warnings;

while (<>) {
    chomp;
    next if /^#/;
    my @F = split /\t/;
    if ($F[2] eq 'LTR_retrotransposon') {
	my $l = $F[4] - $F[3] + 1; 
	my ($id) = ($F[8] =~ /ID=(LTR_retrotransposon\d+)/); 
	my ($sim) = ($F[8] =~ /ltr_similarity=(\d+.\d+)/); 
	print join "\t", $id, $l, $sim, "\n";
    }
}
