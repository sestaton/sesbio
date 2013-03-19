#!/usr/bin/env perl

use v5.12;
use strict;
use warnings;
use autodie qw(open);

my $usage = "$0 annot.tsv\n";
my $infile = shift or die $usage;

my %annot;
my $total_ct = 0;

open(my $in, '<', $infile);

while (<$in>) {
    chomp;
    next if /^Cluster/;
    my @fields = split;
    $total_ct += $fields[1];
    if (scalar @fields == 4) {
	if (exists $annot{$fields[2]}{$fields[3]}) {
	    push @{$annot{$fields[2]}{$fields[3]}}, $fields[1];
	}
	else {
	    $annot{$fields[2]}{$fields[3]} = [$fields[1]];
	}
    }
    else {
	my $fam = $fields[5];
	$fam =~ s/\-\d.*//;
	if (exists $annot{$fields[4]}{$fam}) {
	    push @{$annot{$fields[4]}{$fam}}, $fields[1];
	}
	else {
	    $annot{$fields[4]}{$fam} = [$fields[1]];
	}
    }
}
close($in);

for my $sf (keys %annot) {
    for my $f (reverse sort { scalar @{$annot{$sf}{$a}} <=> scalar @{$annot{$sf}{$b}} } keys %{$annot{$sf}}) {
	my $read_ct = scalar @{$annot{$sf}{$f}};
	my $perc_cov = sprintf("%.12f",$read_ct/$total_ct);
	say join "\t", $sf, $f, $read_ct."/".$total_ct, $perc_cov;
    }
}
