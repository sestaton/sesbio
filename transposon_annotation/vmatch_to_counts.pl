#!/usr/bin/env perl

## This script will take the output from Vmatch 
## and generate a file of counts by repeat type (superfamily)

use 5.010;
use strict;
use warnings;
use autodie;
#use Data::Dump;
use JSON;
use List::Util qw(sum);

my $usage    = "$0 idlist vmatches\n";
my $infile   = shift or die $usage;
my $vmatches = shift or die $usage;
open my $in, '<', $infile;
open my $vm, '<', $vmatches;

my %family_map;

while (my $line = <$in>) {
    chomp $line;
    my ($f, $sf, $source)  = split /\t/, $line;
    next unless defined $sf && defined $f; ## why?
    if ($sf =~ /(\s+)/) {
	$sf =~ s/$1/\_/;
    }
    $f =~ s/\s/\_/;
    if (exists $family_map{$sf}) {
	push @{$family_map{$sf}}, {$f => 0};
    }
    else {
	$family_map{$sf} = [];
    }
}
close $in;

while (my $l = <$vm>) {
    chomp $l;
    $l =~ s/^\s+//;
    my ($ct, $fam_match)  = split /\s/, $l;
    if ($ct > 1) {
	for my $superfamily (keys %family_map) {
	    while (my ($family_index, $family) = each @{$family_map{$superfamily}}) {
		for my $fam_name (keys %$family) {
		    if (exists $family->{$fam_match}) {
			$family_map{$superfamily}[$family_index]{$fam_name} = $ct;
		    }
		}
	    }
	}
    }
}
close $vm;

my %supfamily_match_ct;

for my $supfam (keys %family_map) {
    while (my ($fam_index, $fam) = each @{$family_map{$supfam}}) {
	for my $fam_n (keys %$fam) {
	    if ($fam->{$fam_n} > 0) {
		if (exists $supfamily_match_ct{$supfam}) {
		    push @{$supfamily_match_ct{$supfam}}, $fam->{$fam_n};
		}
		else {
		    $supfamily_match_ct{$supfam} = [];
		}
	    }
	}
    }
    my $sum = sum(@{$supfamily_match_ct{$supfam}});
    if (defined $sum) {
	say join "\t", $supfam, $sum;
    }
}

#dd %supfamily_match_ct;
