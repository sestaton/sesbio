#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;
use autodie;
use Bio::DB::HTS::Faidx;
use Data::Dump::Color;

my $usage = "$0 windows.txt seq";
my $windows = shift or die $usage;
my $seqfile = shift or die $usage;

my $index = Bio::DB::HTS::Faidx->new($seqfile);
open my $w, '<', $windows;

my %bins;
my $binct = 0;
while (my $bline = <$w>) {
    chomp $bline;
    $binct++;
    my ($chr, $start, $end) = split /\t/, $bline;
    #HanXRQChr01 0 1000000 
    #HanXRQChr01 1000000 2000000 
    my ($fct, $btot) = (0, 0);
    my $id = "$chr:$start-$end";
    my ($seq, $len)  = $index->get_sequence(id);
    my $gc_num = $seq =~ tr/GCgc/GCgc/;
    my $gc_content = sprintf(".2f", $gc_num/$len);
    say join "\t", $chr, $start, $end, $gc_content;
}
close $w;
