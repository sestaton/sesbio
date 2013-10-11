#!/usr/bin/env perl

use v5.12;
use strict;
use warnings;
use autodie qw(open);
use Statistics::Descriptive;

my $usage = "$0 infile > outfile\n";
my $infile = shift or die $usage;
my %fams;

open my $in, '<', $infile;
while (<$in>) {
    chomp;
    next if /^Family/;
    my ($fam, $gcov) = split;
    if (exists $fams{$fam}) {
	push @{$fams{$fam}}, $gcov;
    }
    else {
	$fams{$fam} = [$gcov];
    }
}

my $stat = Statistics::Descriptive::Full->new();
my %cv;
say join "\t", "Family", "CV";
for my $ufam (sort { $fams{$a} cmp $fams{$b} } keys %fams) {
    if (scalar @{$fams{$ufam}} > 1) {
	$stat->add_data(@{$fams{$ufam}});
	my $cv = 100*($stat->standard_deviation/$stat->mean);
	#say join "\t", $ufam, $cv;
	$cv{$ufam} = $cv;
    }
    else {
	$cv{$ufam} = 0;
	#say join "\t", $ufam, 0;
    }
}

for my $k (reverse sort { $cv{$a} <=> $cv{$b} } keys %cv) {
    say join "\t", $k, $cv{$k};
}
