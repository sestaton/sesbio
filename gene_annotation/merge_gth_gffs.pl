#!/usr/bin/env perl

## This is for merging GFFs from GenomeThreader run on separate chromosomes
## into a single GFF for training.

use 5.010;
use strict;
use warnings;
use Sort::Naturally;

my @gffs = glob "Ha*gth.gff3";

if (@gffs < 1) {
    say "\nERROR: Expecting GFF files matching 'Ha*gth.gff3'. Exiting.\n";
    exit(1);
}

my ($header, %files);

for my $gff (nsort @gffs) {
    open my $in, '<', $gff;
    while (my $line = <$in>) {
	chomp $line;
	if ($line =~ /^##sequence-region/) {
	    my ($sr, $loc, $s, $e) = split /\s+/, $line;
	    my ($reg) = ($loc =~ /(Ha\d+)/);
	    $header .= join q{ }, $sr, $reg, $s, $e."\n";
	}
    }
    close $in;
}

say "##gff-version 3";
chomp $header;
say $header;

my ($genes, $mrna, $cds) = (0, 0, 0);
for my $gff (nsort @gffs) {
    say STDERR "Working on...$gff";

    open my $in, '<', $gff;
    while (my $line = <$in>) {
	chomp $line;
	next if $line =~ /^#/;
	my @feats = split /\t/, $line;
	my ($reg) = ($feats[0] =~ /(Ha\d+)/);
	if ($feats[2] eq 'gene') {
	    $genes++;
	}
	if ($feats[2] eq 'mRNA') {
	    $mrna++;
	}
	
	$feats[8] =~ s/\s\;\s/\;/g;
	$feats[8] =~ s/\s+$//;
	$feats[8] =~ s/\"//g;
	$feats[8] =~ s/(\;\w+)\s/$1=/g;
	$feats[8] =~ s/\s;/;/;
	$feats[8] =~ s/^(\w+)\s/$1=/;
	$feats[8] =~ s/(?:Target=md5:.*):/Target=$reg:/;
	$feats[8] =~ s/ID=gene\d+/ID=gene$genes/g;
	$feats[8] =~ s/Parent=gene\d+/Parent=gene$genes/g;	
	$feats[8] =~ s/mRNA\d+/mRNA$mrna/g;
	$feats[8] =~ s/CDS\d+/CDS$mrna/g;
	say join "\t", $reg, @feats[1..8];
    }
    close $in;
}
