#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;
use autodie qw(open);
use Bio::DB::HTS::Kseq;
use Data::Dump::Color;

my $usage = "$0 genes.fas blast.bln 1> filtered_list.txt 2> filtered_results.txt";
my $fasta = shift or die $usage;
my $blast = shift or die $usage;
my $frac  = 0.80; # fraction coverage
my $pid   = 90;   # percent identity threshold

my %genes;
my $kseq = Bio::DB::HTS::Kseq->new($fasta);
my $iter = $kseq->iterator;
while (my $seqobj = $iter->next_seq) {
    my $id = $seqobj->name;
    my $seq = $seqobj->seq;
    $genes{$id} = length($seq);
}
my $genect = keys %genes;
#dd \%genes and exit;
#say join "\t", "#Query", "Hit", "Percent_ID", "HSP_len", "Num_mismatch", "Num_gaps", 
#                          "Query_start", "Query_end", "Hit_start", "Hit_end", "E-value", "Bit_score";

my %removed;
open my $bin, '<', $blast;
while (my $line = <$bin>) {
    chomp $line;
    my @f = split /\t/, $line;
    if (exists $genes{$f[0]}) {
	if ($f[3]/$genes{$f[0]} >= $frac) {
	    $removed{$f[0]} = sprintf("%.2f",$f[3]/$genes{$f[0]})
		unless exists $removed{$f[0]};
	    #say join q{ }, $f[0], $genes{$f[0]}, $f[3], sprintf("%.2f",$f[3]/$genes{$f[0]});
	}
    }
    else {
	warn $f[0]," not found";
    }
}
close $bin;

my $rmct = keys %removed;
say STDERR "$genect total genes";
say STDERR "$rmct genes over threshold";
say STDERR ($genect-$rmct)," filtered genes";

for my $gene (keys %genes) {
    unless (exists $removed{$gene}) {
	say $gene;
    }
}
