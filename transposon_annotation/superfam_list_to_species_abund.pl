#!/usr/bin/env perl

## NB: This script will grab all the Transposome annotation summary files in the working directory
##     and calculate the transposon abundance for each species (where each file is assumed to 
##     be from a different species) and write the results in tab-delimited format to STDOUT.
##
##     The number of expected species is hard-coded at 15 on line 21. This is a test to make sure you are
##     not grabbing more files than expected. Change as need for the number of species.

use 5.010;
use strict;
use warnings;
use autodie qw(open);
use Statistics::Descriptive;
#use Data::Dump::Color;

my $usage = "perl $0\n";

my @files = glob "*summary.tsv";

die "\nERROR: Could not find and summary files." unless scalar @files == 15; ## 15 species for me, change as needed

my %typemap = (
    'Academ'       => { transposable_element => 'dna_transposon'          },
    'Ambal'        => { transposable_element => 'non-ltr_retrotransposon' },
    'BEL'          => { transposable_element => 'ltr_retrotransposon'     },
    'CRE'          => { transposable_element => 'dna_transposon'          },
    'Copia'        => { transposable_element => 'ltr_retrotransposon'     },
    'Crack'        => { transposable_element => 'dna_transposon'          },
    'Crypton'      => { transposable_element => 'dna_transposon'          },
    'DIRS'         => { transposable_element => 'ltr_retrotransposon'     },
    'Daphne'       => { transposable_element => 'non-ltr_retrotransposon' },
    'ERV1'         => { transposable_element => 'endogenous_retrovirus'   },
    'ERV2'         => { transposable_element => 'endogenous_retrovirus'   },
    'EnSpm'        => { transposable_element => 'dna_transposon'          },
    'Ginger2/TDD'  => { transposable_element => 'dna_transposon'          },
    'Gypsy'        => { transposable_element => 'ltr_retrotransposon'     },
    'Harbinger'    => { transposable_element => 'dna_transposon'          },
    'Helitron'     => { transposable_element => 'dna_transposon'          },
    'ISL2EU'       => { transposable_element => 'dna_transposon'          },
    'Jockey'       => { transposable_element => 'non-ltr_retrotransposon' },
    'Kiri'         => { transposable_element => 'non-ltr_retrotransposon' },
    'L1'           => { transposable_element => 'non-ltr_retrotransposon' },
    'L2'           => { transposable_element => 'non-ltr_retrotransposon' },
    'Mariner/Tc1'  => { transposable_element => 'dna_transposon'          },
    'MuDR'         => { transposable_element => 'dna_transposon'          },
    'NeSL'         => { transposable_element => 'non-ltr_retrotransposon' },
    'Penelope'     => { transposable_element => 'non-ltr_retrotransposon' },
    'Polinton'     => { transposable_element => 'dna_transposon'          },
    'Proto1'       => { transposable_element => 'non-ltr_retrotransposon' },
    'Proto2'       => { transposable_element => 'non-ltr_retrotransposon'},
    'R1'           => { transposable_element => 'non-ltr_retrotransposon' },
    'R2'           => { transposable_element => 'non-ltr_retrotransposon' },
    'RTE'          => { transposable_element => 'non-ltr_retrotransposon' },
    'RTEX'         => { transposable_element => 'non-ltr_retrotransposon' },
    'SINE2/tRNA'   => { transposable_element => 'non-ltr_retrotransposon' },
    'Sola'         => { transposable_element => 'dna_transposon' },
    'Tad1'         => { transposable_element => 'non-ltr_retrotransposon' },
    'Transib'      => { transposable_element => 'dna_transposon' },
    'Tx1'          => { transposable_element => 'non-ltr_retrotransposon' },
    'hAT'          => { transposable_element => 'dna_transposon' },
    );

my %species_types;
my %classes;

for my $file (@files) {
    my ($sp, $rest) = split /\_/, $file, 2;
    open my $in, '<', $file;
    while (<$in>) {
	chomp;
	#Superfamily Family Mean SD
	next if /^Superfamily/;
	my @f = split;
	if (exists $typemap{$f[0]}) {
	    my ($type, $class) = map { $_, $typemap{$f[0]}->{$_} } keys %{$typemap{$f[0]}};    
	    push @{$classes{$class}}, $f[2];
	    push @{$species_types{$sp}{$f[0]}}, { $f[1] => $f[2] };
	}
	else {
	    say "Superfamily NOT defined =====> ", $f[0];
	}
    }
    close $in;
}

#for my $k (sort keys %classes) {
#    my $stat = Statistics::Descriptive::Full->new;
#    $stat->add_data(@{$classes{$k}});
#    my $mean = $stat->mean;
#    my $sd   = $stat->standard_deviation;
#    say join "\t", $k, $mean, $sd;
#}

#dd \%species_types;
for my $species (keys %species_types) {
    for my $sfam (keys %{$species_types{$species}}) {
	if (exists $typemap{$sfam} ) {
	    my ($type, $class ) = map { $_, $typemap{$sfam}->{$_} } keys %{$typemap{$sfam}};
	    for my $fam (@{$species_types{$species}{$sfam}}) {
		for my $famname (keys %$fam) {
		    say join "\t", $species, $type, $class, $sfam, $famname, $fam->{$famname};
		}
	    }
	}
    }
}

