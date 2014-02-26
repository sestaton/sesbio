#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;

my %res;

while (<>) {
    chomp;
    my @f = split;
    if (exists $res{$f[0]}) {
	push @{$res{$f[0]}}, { $f[3] => $f[5] };
    }
    else {
	$res{$f[0]} = [ {$f[3] => $f[5] } ];
    }
}

my %sfamtot;

for my $species (keys %res) {
    for my $sfamh (@{$res{$species}}) {
	my $sfam_tot = 0;
	for my $sfam (keys %$sfamh) {
	   $sfamtot{$sfam} += $sfamh->{$sfam};
	}
    }
    for my $sfamname (reverse sort { $sfamtot{$a} <=> $sfamtot{$b} } keys %sfamtot) {
	say join "\t", $species, $sfamname, $sfamtot{$sfamname};
    }
    %sfamtot = ();
}
	
