#!/usr/bin/env perl


## TODO: Change names for printing out full Genus name in plot

use 5.014;
use utf8;
use strict;
use warnings;
use warnings FATAL => "utf8";
use autodie qw(open);
use Data::Dump qw(dd);
use charnames qw(:full :short);

my $usage = "perl $0 outfile\n";
my $outfile = shift or die $usage;

my @files = glob("*summary.tsv");

die "\nERROR: Could not get files: $!"
    unless scalar @files > 0;

my %df;
my %fams;
my %sph;

my $cvalh = {
    Phoebanthus => 4267295897,
    Ageratina   => 1746269472,
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
    Senecio     => 2045909989,
    TKS         => 2582325776,
};

for my $file (@files) {
    my ($species) = split(/\_/, $file, 2);
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

    open my $in, '<', $file;
    while (<$in>) {
	chomp;
	next if /^ReadNum/;
	my @f = split "\t";
	$fams{$f[2]} = 1;
	my $key = mk_key($f[2], $species);
	$sph{$key} = 1;
	$df{$species}{$f[2]} = $f[5];
    }
    close $in;
}

open my $out, '>', $outfile;
my $fam_ct = scalar keys %fams;
print $out "Species\t", join "\t", (sort keys %fams), "\n";
keys %fams;

my %species_map = ( 
    Ager      => 'Ageratina',
    Ann1238   => 'Hannuus',
    Harg      => 'Hargophyllus',
    Hport     => 'Hporterii',
    Hteph     => 'Hniveus',
    Hvert     => 'Hverticillatus',
    CP        => 'Centrapallus',
    Calyc     => 'Nasanthus',
    Dasy      => 'Fulcaldea',
    Gerb      => 'Gerbera',
    Gnaph     => 'Pseudognaphalium',
    Saff      => 'Carthamus',
    Sene      => 'Senecio',
    TKS       => 'Taraxacum',
    Phoeb     => 'Phoebanthus',
    );

##TODO: Fix bug with creating format that can't be read in R without reformatting in excel
for my $species (sort keys %df) {
    print $out $species_map{$species};
    for my $f (sort keys %fams) {
	my $key = mk_key($f, $species);
	if (exists $sph{$key} && exists $df{$species}{$f} && exists $cvalh->{$species}) {
	    my $bpsize = $df{$species}{$f} * $cvalh->{$species};
	    my $rounded = int($bpsize + $bpsize/abs($bpsize*2));
	    print $out "\t$rounded";
	}
	else {
	    print $out "\t0";
	}
    }
    print $out "\n";
}
close $out;

## subs
sub mk_key { join "\N{INVISIBLE SEPARATOR}", map { $_ // " " } @_ }

sub mk_vec { split "\N{INVISIBLE SEPARATOR}", shift }
