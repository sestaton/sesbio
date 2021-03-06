#!/usr/bin/env perl

## This script is for calculating statistical differences for certain TE types between species.
## The assumption is that there are annotation summary files from Transposome in the working
## directory that would correspond to different species.

use 5.010;
use strict;
use warnings;
use autodie qw(open);
use Statistics::Descriptive;
use Data::Dump qw(dd);

my $usage = "perl $0\n";

my @files = glob "*summary.tsv";

die "\nERROR: Could not find all summary files." 
    unless scalar @files == 15; ## corresponds to 15 species, change as needed but use as a sanity check

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

my %allspecies;
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
	    my ($type, $class ) = map { $_, $typemap{$f[0]}->{$_} } keys %{$typemap{$f[0]}};    
	    push @{$classes{$sp}{$class}}, $f[2];
	}
	else {
	    say "Time for more work: DEBUG! ====> ", $f[0];
	}
    }
    close $in;
}

#dd \%classes and exit;

my %allsp;
for my $species (sort keys %classes) {
    for my $cls (keys %{$classes{$species}}) {
	my $stat = Statistics::Descriptive::Full->new;

	if (scalar @{$classes{$species}{$cls}} == 1) {
	    push @{$classes{$species}{$cls}}, (0, 0);
	}
	elsif (scalar @{$classes{$species}{$cls}} == 2) {
	    push @{$classes{$species}{$cls}}, 0;
	}

	$stat->add_data(@{$classes{$species}{$cls}});
	my $mean = $stat->mean;
	my $sd   = $stat->standard_deviation;
	my $sum  = $stat->sum;
	say join "\t", $species, $cls, $sum, $mean, $sd;
    }
    #undef $stat;
}
