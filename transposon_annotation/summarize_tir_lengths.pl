#!/usr/bin/env perl

use 5.020;
use strict;
use warnings;
use autodie;
use Statistics::Descriptive;
use Sort::Naturally;
use Bio::Tools::GFF;
use List::Util qw(any);
use Data::Dump;

my $usage = "$0 gff";
my $gff = shift or die $usage;

open my $in, '<', $gff;
while (<$in>) {
    chomp;
    if (/^#/) {
	say $_;
    }
    else {
	last;
    }
}
close $in;

my $ltrrt = 0;
my @rt = qw(rve rvt rvp gag chromo);
my $gffio = Bio::Tools::GFF->new( -file => $gff, -gff_version => 3 );

my ($start, $end, $region, %features);
while (my $feature = $gffio->next_feature()) {
    if ($feature->primary_tag eq 'repeat_region') {
	my @string = split /\t/, $feature->gff_string;
	($region) = ($string[8] =~ /ID=?\s+?(repeat_region\d+)/);
	($start, $end) = ($feature->start, $feature->end);
    }
    next $feature unless defined $start && defined $end;
    if ($feature->primary_tag ne 'repeat_region') {
	if ($feature->start >= $start && $feature->end <= $end) {
	    push @{$features{$region.".".$start.".".$end}}, join "||", split /\t/, $feature->gff_string;
	}
    }
}

#dd \%feature;
my ($tirct, $filtct) = (0, 0);
for my $tir (nsort keys %features) {
    $tirct++;
    for my $feat (@{$features{$tir}}) {
	my @feats = split /\|\|/, $feat;
	if ($feats[2] eq 'protein_match') {
	    my ($type, $pdom) = ($feats[8] =~ /(name) ("?\w+"?)/);
	    $pdom =~ s/"//g;
	    my $dom = lc($pdom);
	    if ($dom =~ /rve|rvt|rvp|gag|chromo|rnase|athila|zf/) {
		delete $features{$tir};
	    }
	}
    }
}

my @lengths;

for my $tir (keys %features) {
    my ($rreg, $s, $e) = split /\./, $tir;
    my $len = ($e - $s);
    push @lengths, $len;
    my $region = @{$features{$tir}}[0];
    my ($loc, $source) = (split /\|\|/, $region)[0..1];
    say join "\t", $loc, $source, 'repeat_region', $s, $e, '.', '?', '.', "ID=$rreg";
    for my $feat (@{$features{$tir}}) {
        my @feats = split /\|\|/, $feat;
	$feats[8] =~ s/\s\;\s/\;/g;
	$feats[8] =~ s/\s+/=/g;
	say join "\t", @feats;
    }
}

my $stat = Statistics::Descriptive::Full->new;
$stat->add_data(@lengths);
my $min  = $stat->min;
my $max  = $stat->max;
my $mean = $stat->mean;

say STDERR join q{ }, $min, $max, $mean;
#$filtct = (keys %feature);
#say STDERR $tirct;
#say STDERR $filtct;
