#!/usr/bin/env perl

use v5.10;
use strict;
use warnings;

my %pair;

while (<>) {
    chomp;
    my @f = split;
    next if $f[0] eq $f[1];
    my @i = split /\_/, $f[0];
    my @j = split /\_/, $f[1];
    #say "$i[0] => $j[0]";
    unless ($i[0] eq $j[0]) {
	my $p = join "||", $f[0], $f[1];
	unless (exists $pair{$p}) {
	    if ($f[3] > 1000) {
		say join "\t", @f;
	    }
	}
	$pair{$p} = 1;
    }
}
