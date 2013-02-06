#!/usr/bin/env perl

use strict;
use warnings;
use autodie qw(open);
use lib qw(/home/jmblab/statonse/apps/perlmod/Data-Dump-1.21/blib/lib);
use Data::Dump qw(dd);

#my $usage = "$0 infile outfile\n";
#my $infile = shift or die $usage;
#my $outfile = shift or die $usage;

#open(my $in, '<', $infile);

my $matches = {};

$matches->{'transposable_element'}{'dna_transposon'} = ['Mariner/Tc1', 'hAT', 'MuDR', 'EnSpm', 'piggyBac', 'P', 'Merlin', 'Harbinger', 'Transib', 'Novosib', 'Helitron', 'Polinton', 'Kolobok', 'ISL2EU', 'Crypton', 'Sola', 'Zator', 'Ginger/1', 'Ginger2/TDD', 'Academ', 'Zisupton', 'IS3EU'];

$matches->{'transposable_element'}{'ltr_retrotransposon'} = ['Gypsy', 'Copia', 'BEL', 'DIRS'];

$matches->{'transposable_element'}{'endogenous_retrovirus'} = ['ERV1', 'ERV2', 'ERV3', 'Lentivirus', 'ERV4'];

$matches->{'transposable_element'}{'non-ltr_retrotransposon'} = [{'SINE' => ['SINE1/7SL','SINE2/tRNA','SINE3/5S','SINE4']},'CRE', 'NeSL', 'R4', 'R2', 'L1', 'RTE', 'I', 'Jockey', 'CR1', 'Rex1', 'RandI', 'Penelope', 'Tx1', 'RTEX', 'Crack', 'Nimb', 'Proto1', 'Proto2', 'RTETP', 'Hero', 'L2', 'Tad1', 'Loa', 'Ingi', 'Outcast', 'R1', 'Daphne', 'L2A', 'L2B', 'Ambal', 'Vingi', 'Kiri'];

$matches->{'simple_repeat'}{'Satellite'} = ['SAT', 'MSAT'];

$matches->{'pseudogene'} = ['rRNA', 'tRNA', 'snRNA'];

$matches->{'integrated_virus'} = ['DNA_Virus', 'Caulimoviridae'];

dd $matches;
