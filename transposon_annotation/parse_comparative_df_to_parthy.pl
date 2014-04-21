#!/usr/bin/env perl

#NB: This is for staton et al. 2014
use 5.014;
use strict;
use warnings;
use Statistics::Descriptive;
use autodie qw(open);
use Data::Dump qw(dd);

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
    my @f = split '\t';
    my $species = shift @f;
    $species =~ s/\s/\_/g;
    $fams{$species} = \@f;
}

for my $sp (keys %fams) {
    next if $sp eq '';
    if (exists $cvalh->{$sp}) {
	my @sp_abund = grep { $_ } @{$fams{$sp}};
	my $outfile = $sp."_parthy_in.txt";
	open my $out, '>', $outfile;
	say $out join "\n", @sp_abund;
	close $out;
    }
    else {
	say "$sp not found in map";
    }
}

