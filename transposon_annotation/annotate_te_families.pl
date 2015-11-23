#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;
use autodie;
use Data::Dump::Color;
use List::Util qw(max);

my $usage = "$0 table blast";
my $table = shift or die $usage;
my $blast = shift or die $usage;

my (%matches, %te_cts);
open my $tin, '<', $table;
while (<$tin>) {
    chomp $tin;
    next if /^Family/;
    my @f = split /\t/;
    $f[2] =~ s/\s+$//;
    $f[1] =~ s/.*\///;
    $te_cts{$f[1]}{$f[0]} = $f[2];
}
close $tin;

open my $in, '<', $blast;
while (<$in>) {
    chomp;
    my @f = split /\t/;
    $f[0] =~ s/_Chr.*//;
    if ($f[0] =~ /RLC/) {
	push @{$matches{Copia}{$f[0]}}, $f[1];
    }
    elsif ($f[0] =~ /RLG/) {
        push @{$matches{Gypsy}{$f[0]}}, $f[1];
    }
    elsif ($f[0] =~ /RLX/) {
	push @{$matches{Unknown}{$f[0]}}, $f[1];
    }
}
close $in;

for my $sf (keys %matches) {
    for my $te (keys %{$matches{$sf}}) {
	my %seen;
	my $te_count = map { $seen{$_}++ } @{$matches{$sf}{$te}};
	my $max = max(values %seen);
	my %rev = reverse %seen;
	my $k   = $rev{$max};
	$k =~ s/_I$|_LTR$//;
	$k =~ s/I$|LTR$//;
	if (exists $te_cts{$sf}{$k}) {
	    say join q{ }, $sf, $te, scalar(@{$matches{$sf}{$te}}), $k, $te_cts{$sf}{$k};
	}
	else {
	    warn "\nERROR: $sf -> $k not found in map";
	}
    }
}
#dd \%matches;
