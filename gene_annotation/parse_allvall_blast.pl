#!/usr/bin/env perl

# find similarity in a genomic assembly

use 5.010;
use strict;
use warnings;

my %pair;

while (<>) {
    chomp;
    my @f = split;
    next if $f[0] eq $f[1];    # skip self hits
    my @i = split /\_/, $f[0];
    my @j = split /\_/, $f[1];
    #say "$i[0] => $j[0]";
    unless ($i[0] eq $j[0]) {
	my $p = join "||", $f[0], $f[1];
	unless (exists $pair{$p}) { # only keep unique matches
	    if ($f[3] > 1000) {     # find alignments over 1kb
		say join "\t", @f;
	    }
	}
	$pair{$p} = 1;
    }
}
