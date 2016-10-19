#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;
use autodie;
use List::Util qw(sum);
use List::MoreUtils qw(uniq);
use Statistics::Descriptive;
use Data::Dump::Color;

my $usage = "$0 filtered.tsv";
my $match = shift or die $usage;

my %h;
say join "\t", 'Element_with_gene', 'AlignedCount', 'UniqProteins';

open my $m, '<', $match;
while (my $line = <$m>) {
    chomp $line;
    my @f = split /\t/, $line;
    #Type Element_with_gene Protein PID AlignLength QueryLength AlignPercent Description
    #Gypsy RLG_singleton_family222_LTR_retrotransposon312_HanXRQChr00c0125_1667_12130 HanXRQCPg0579741 99.26 408 500 81.60 Putative..
    push @{$h{$f[1]}}, $f[2];
}
close $m;

for my $k (reverse sort { @{$h{$a}} <=> @{$h{$b}} } keys %h) {
    #my $stat = Statistics::Descriptive::Full->new;
    #$stat->add_data(@{$h{$k}});
    #my $ct = $stat->count;
    my $ct = @{$h{$k}};
    my $uniq = uniq(@{$h{$k}});
    say join "\t", $k, $ct, $uniq;
}
