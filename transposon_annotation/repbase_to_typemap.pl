#!/usr/bin/env perl

use 5.014; ## This is important, we need at least v5.12 for this script
use strict;
use warnings;
use autodie qw(open);
use Data::Dump qw(dd);
use JSON;

## NB: The input is a RepBase ID list (grep ">" repbase.fasta | sed 's/>//g' > idlist).
##     Though, there is code in transposome to perform this same task without creating the ID list manually.

my $usage = "$0 idlist\n";
my $infile = shift or die $usage;
open my $in, '<', $infile;

my $matches = build_repbase_hash();

my %family_map;

while (my $line = <$in>) {
    chomp $line;
    my ($f, $sf, $source)  = split /\t/, $line;
    next unless defined $sf && defined $f; ## why?
    if ($sf =~ /(\s+)/) {
	$sf =~ s/$1/\_/;
    }
    $f =~ s/\s/\_/;
    if (exists $family_map{$sf}) {
	push @{$family_map{$sf}}, $f;
    }
    else {
	$family_map{$sf} = [];
    }
}
close $in;

for my $type (keys %$matches) {
    unless ($type eq 'pseudogene' || $type eq 'integrated_virus') {
	for my $class (keys %{$matches->{$type}}) {
	    while ( my ($superfam_index, $superfam) = each @{$matches->{$type}{$class}} ) {
		for my $superfam_h (keys %$superfam) {
		    my $superfam_cp = lc($superfam_h);
		    for my $mapped_fam (keys %family_map) {
			my $mapped_fam_cp = lc($mapped_fam);
			if (length($superfam_cp) > 1 && length($mapped_fam_cp) > 1) {
			    if ($mapped_fam_cp =~ /sine/i && $superfam_cp =~ /sine/i) {
				while (my ($sine_fam_index, $sine_fam_h) = each @{$superfam->{$superfam_h}}) {
				    for my $sine_fam_mem (keys %$sine_fam_h) {
					if ($sine_fam_mem =~ /$mapped_fam_cp/i && $mapped_fam_cp =~ /^(?!sine$)/) {
					    push @{$matches->{$type}{$class}[$superfam_index]{$superfam_h}[$sine_fam_index]{$sine_fam_mem}}, 
					    $family_map{$mapped_fam};
					}
				    }
				}
			    } 
			    elsif ($mapped_fam_cp =~ /$superfam_cp/) {
				push @{$matches->{$type}{$class}[$superfam_index]{$superfam_h}}, $family_map{$mapped_fam};
			    }
			    elsif (length($mapped_fam_cp) == 1 && length($superfam_cp) == 1) {
				if ($mapped_fam_cp =~ /$superfam_cp/) {
				    push @{$matches->{$type}{$class}[$superfam_index]{$superfam_h}}, $family_map{$mapped_fam};
				}
			    }
			}
		    }
		}
	    } 
	} 
    }
}

#dd $matches;
my $json = JSON->new->utf8->space_after->encode($matches);

my $hash = JSON->new->utf8->space_after->decode($json);
#dd $hash;

#
# subs
#
sub build_repbase_hash {
    my $matches = {};
    
    $matches->{'transposable_element'}{'dna_transposon'} = [{'Mariner/Tc1' => []}, {'hAT' => []}, 
							    {'MuDR' => []}, {'EnSpm' => []}, 
							    {'piggyBac' => []}, {'P' => []}, 
							    {'Merlin' => []}, {'Harbinger' => []}, 
							    {'Transib' => []}, {'Novosib' => []}, 
							    {'Helitron' => []}, {'Polinton' => []}, 
							    {'Kolobok' => []}, {'ISL2EU' => []}, 
							    {'Crypton' => []}, {'Sola' => []}, 
							    {'Zator' => []}, {'Ginger/1' => []}, 
							    {'Ginger2/TDD' => []}, {'Academ' => []}, 
							    {'Zisupton' => []}, {'IS3EU' => []}];
    
    $matches->{'transposable_element'}{'ltr_retrotransposon'} = [{'Gypsy' => []}, {'Copia' => []}, 
								 {'BEL' => []}, {'DIRS' => []}];
    
    $matches->{'transposable_element'}{'endogenous_retrovirus'} = [{'ERV1' => []}, {'ERV2' => []}, 
								   {'ERV3' => []}, {'Lentivirus' => []}, 
								   {'ERV4' => []}];
    
    $matches->{'transposable_element'}{'non-ltr_retrotransposon'} = [{'SINE' => [{'SINE1/7SL' => []}, {'SINE2/tRNA' => []},
										 {'SINE3/5S' => []},{'SINE4' => []}]},
								     {'CRE' => []}, {'NeSL' => []}, 
								     {'R4' => []}, {'R2' => []}, 
								     {'L1' => []}, {'RTE' => []}, 
								     {'I' => []}, {'Jockey' => []}, 
								     {'CR1' => []}, {'Rex1' => []}, 
								     {'RandI' => []}, {'Penelope' => []}, 
								     {'Tx1' => []}, {'RTEX' => []}, 
								     {'Crack' => []}, {'Nimb' => []}, 
								     {'Proto1' => []}, {'Proto2' => []}, 
								     {'RTETP' => []}, {'Hero' => []}, 
								     {'L2' => []}, {'Tad1' => []}, 
								     {'Loa' => []}, {'Ingi' => []}, 
								     {'Outcast' => []}, {'R1' => []}, 
								     {'Daphne' => []}, {'L2A' => []}, 
								     {'L2B' => []}, {'Ambal' => []}, 
								     {'Vingi' => []}, {'Kiri' => []}];

    $matches->{'simple_repeat'}{'Satellite'} = [{'SAT' => []}, {'MSAT' => []}];

    $matches->{'pseudogene'} = [{'rRNA' => []}, {'tRNA' => []}, {'snRNA' => []}];

    $matches->{'integrated_virus'} = [{'DNA_Virus' => []}, {'Caulimoviridae' => []}];
    
    return $matches;
}
