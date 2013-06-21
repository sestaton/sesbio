#!/usr/bin/env perl


## TODO: Change names for printing out full Genus name in plot

use v5.14;
use utf8;
use strict;
use warnings;
use warnings FATAL => "utf8";
use autodie qw(open);
use Data::Dump qw(dd);
use charnames qw(:full :short);

my $usage = "perl $0 outfile\n";
my $outfile = shift or die $usage;

my @files = glob("*summary_edit.tsv");

die "\nERROR: Could not get files: $!"
    unless scalar @files > 0;

my %df;
my %fams;
my %sph;

for my $file (@files) {
    my ($species) = split(/\_/, $file, 2);
    open my $in, '<', $file;
    while (<$in>) {
	chomp;
	next if /^ReadNum/;
	my @f = split "\t";
	if (exists $df{$species}) {
	    push @{$df{$species}}, $f[2];
	}
	else {
	    $df{$species} = [ $f[2] ];
	}
    }
    close $in;
}

my %species_map = ( 
    Ageratina => 'Ageratina',
    Ann1238   => 'Helianthus',
    CP        => 'Centrapallus',
    Calyc     => 'Nasanthus',
    Dasy      => 'Fulcaldea',
    Gerb      => 'Gerbera',
    Gnaph     => 'Pseudognaphalium',
    Saff      => 'Carthamus',
    Senecio   => 'Senecio',
    TKS       => 'Taraxacum',
    );


say join "\t", "Species", "Families";
for my $sp (sort keys %df) {
    my $sp_fam_ct = scalar @{$df{$sp}};
    say join "\t", $species_map{$sp}, $sp_fam_ct;
}
