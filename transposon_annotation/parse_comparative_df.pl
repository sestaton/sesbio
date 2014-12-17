#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;
use Statistics::Descriptive;
use autodie    qw(open);
#use Data::Dump qw(dd);

my $cvalh = {
    'Phoebanthus_tenuifolius'          => 4267295897,
    'Conoclinium_coelestinum'          => 1746269472,
    'Helianthus_annuus'                => 3384161947,
    'Helianthus_argophyllus'           => 4174346891,
    'Helianthus_porteri'               => 4330738740,
    'Helianthus_niveus_ssp._tephrodes' => 4192677026,
    'Helianthus_verticillatus'         => 2278002736,
    'Centrapalus_pauciflorus'          => 3125365235,
    'Nasanthus_patagonicus'            => 3962340892,
    'Fulcaldea_stuessyi'               => 4182557218,
    'Gerbera_hybrida'                  => 3861919879,
    'Pseudognaphalium_helleri'         => 2920317131,
    'Carthamus_tinctorius'             => 2405291468,
    'Senecio_vulgaris'                 => 2045909989,
    'Taraxacum_kok-saghyz'             => 2582325776,
};

my $famct   = 0;
my $famsize = 0;
my %fams;

while (<>) {
    chomp;
    next if /^\s+|^\t/;
    my @f = split '\t';
    my $species = shift @f;
    $species =~ s/\s/\_/g;
    $fams{$species} = \@f;
}

my @all_stats;
my $statsall = Statistics::Descriptive::Full->new;
for my $sp (keys %fams) {
    if (exists $cvalh->{$sp}) {
	my @fam_perc = map  { ($_ / $cvalh->{$sp}) * 100  } 
	               grep { $_ } @{$fams{$sp}};

	my $stat = Statistics::Descriptive::Full->new;
	$statsall->add_data(@fam_perc);
	$stat->add_data(@fam_perc);
	my $mean = $stat->mean;
	my $famnum = @fam_perc;
	push @all_stats, $famnum;
	say join "\t", $sp, $mean, $famnum;
	undef $stat;
    }
    else {
	say "$sp not found in map";
    }
}

my $grandmean = $statsall->mean;
my $sizestats = Statistics::Descriptive::Full->new;
$sizestats->add_data(@all_stats);
my $famnummean = $sizestats->mean;
say '';
say join "\t", "FamSizeMean", "FamNumMean";
say join "\t", $grandmean, $famnummean;
