#!/usr/bin/env perl

use 5.014;
use strict;
use warnings;
use autodie qw(open);
use Data::Dump qw(dd);
use Statistics::Descriptive;

my @summaries = glob "*summary.tsv";
if (scalar @summaries < 1) {
    die "ERROR: No files found!\n";
}

my %sfam_fam;

my $cvalh = {
    Phoeb       => 4267295897,
    Ager        => 1746269472,
    Ann1238     => 3384161947,
    Harg        => 4174346891,
    Hport       => 4330738740,
    Hteph       => 4192677026,
    Hvert       => 2278002736,
    CP          => 3125365235,
    Calyc       => 3962340892,
    Dasy        => 4182557218,
    Gerb        => 3861919879,
    Gnaph       => 2920317131,
    Saff        => 2405291468,
    Sene        => 2045909989,
    TKS         => 2582325776,
};

for my $file (@summaries) {
    my ($species, $rest) = split /\_/, $file, 2; 
    $species = 'Ager'        if $species =~ /ager/;
    $species = 'Ann1238'     if $species =~ /ann1238/;
    $species = 'CP'          if $species =~ /cp/;
    $species = 'Calyc'       if $species =~ /calyc/;
    $species = 'Dasy'        if $species =~ /dasy/;
    $species = 'Gerb'        if $species =~ /gerb/;
    $species = 'Gnaph'       if $species =~ /gnaph/;
    $species = 'Harg'        if $species =~ /harg/;
    $species = 'Hport'       if $species =~ /hport/;
    $species = 'Hteph'       if $species =~ /hteph/;
    $species = 'Hvert'       if $species =~ /hvert/;
    $species = 'Phoeb'       if $species =~ /phoeb/;
    $species = 'Saff'        if $species =~ /saff/;
    $species = 'Sene'        if $species =~ /sene/;
    $species = 'TKS'         if $species =~ /tks/;

    open my $f, '<', $file;
    while (<$f>) {
	chomp;
	next if /^ReadNum/;
	my @res = split;
	$sfam_fam{$species}{$res[1]}{$res[2]} = $res[5];
    }
}

my @fam_size;

my %sp_sum;
my $famct = 0;
#dd \%sfam_fam; exit;
say join "\t", "Species", "Families", "Family_mean", "Genome_size";

for my $spec (keys %sfam_fam) {
    my $stat = Statistics::Descriptive::Full->new;
    for my $sf (keys %{$sfam_fam{$spec}}) { 
	for my $fam (keys %{$sfam_fam{$spec}{$sf}}) {
	    $famct++;
	    $stat->add_data($sfam_fam{$spec}{$sf}{$fam});
	}
	my $mean = $stat->mean;
    }
    my $mean = $stat->mean;
    say join "\t", $spec, $famct, sprintf("%.2f",$mean*100), sprintf("%.0f",$cvalh->{$spec}/1000000);
    $famct = 0;
}
