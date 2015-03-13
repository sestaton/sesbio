#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;
use autodie;
use Data::Dump;
use JSON;
use List::MoreUtils qw(first_index);

my $usage = "$0 fasta\n";
my $infile = shift or die $usage;
open my $in, '<', $infile;

my $matches = build_repbase_hash();
my $repeats = map_repeat_types($matches);
#dd $map and exit;

my %family_map;

while (my $line = <$in>) {
    chomp $line;
    if ($line =~ /^>/) {
	$line =~ s/>//;
	my ($f, $sf, $source)  = split /\t/, $line;
	next unless defined $sf && defined $f; ## why?
	if ($sf =~ /(\s+)/) {
	    $sf =~ s/$1/\_/;
	}
	$f =~ s/\s/\_/;
	push @{$family_map{$sf}}, $f;
    }
}
close $in;

#dd \%family_map and exit;
#my $map;
for my $mapped_sfam (keys %family_map) {
    my $mapped_sfam_cp = lc($mapped_sfam);
    for my $mapped_fam (@{$family_map{$mapped_sfam}}) {
	#$map = map_all_terms($repeats, $mapped_sfam_cp, $mapped_fam);
	for my $class (keys %$repeats) {
	    for my $sfamh (@{$repeats->{$class}}) {
		my $sfam_index = first_index { $_ eq $sfamh } @{$repeats->{$class}};
		for my $sfamname (keys %$sfamh) {
		    if (lc($sfamname) eq $mapped_sfam_cp) {
			push @{$repeats->{$class}[$sfam_index]{$sfamname}}, $mapped_fam;
		    }
		}
	    }
	}

    }
}
#exit;

#dd $repeats and exit;
#my $json = JSON->new->utf8->space_after->encode($matches);

#my $hash = JSON->new->utf8->space_after->decode($json);
#dd $hash;

#
# subs
#
sub map_all_terms {
    my ($map, $superfamily, $family) = @_;
    for my $class (keys %$map) {
	for my $sfamh (@{$map->{$class}}) {
	    my $sfam_index = first_index { $_ eq $sfamh } @{$map->{$class}};
	    for my $sfamname (keys %$sfamh) {
		if (lc($sfamname) eq $superfamily) {
		    push @{$map->{$class}[$sfam_index]{$sfamname}}, $family;
		}
	    }
	}
    }
    return $map;
}

sub map_repeat_types {
    my ($matches) = @_;

    my %repeats;

    for my $type (keys %$matches) { 
	if ($type eq 'transposable_element') { 
	    for my $tes (keys %{$matches->{$type}}) {
		if ($tes eq 'dna_transposon') {
		    $repeats{'dna_transposon'} = $matches->{$type}{$tes};
		}
                elsif ($tes eq 'ltr_retrotransposon') {
                    $repeats{'ltr_retrotransposon'} = $matches->{$type}{$tes};
                }
                elsif ($tes eq 'non-ltr_retrotransposon') {
		    $repeats{'non-ltr_retrotransposons'} = $matches->{$type}{$tes};
                }
                elsif ($tes eq 'endogenous_retrovirus') {
                    $repeats{'endogenous_retrovirus'} = $matches->{$type}{$tes};
                }
            }
	}        
	elsif ($type eq 'simple_repeat') { 
            for my $subtype (keys %{$matches->{$type}}) {
                if ($subtype eq 'Satellite') {
                    $repeats{'satellite'} = $matches->{$type}{$subtype};
                }
            }
        }
        elsif ($type eq 'pseudogene') { 
            $repeats{'pseudogene'} = $matches->{$type}
        }
        elsif ($type eq 'integrated_virus') { 
            $repeats{'integrated_virus'} = $matches->{$type}
        }
    }
    return \%repeats;
}

sub build_repbase_hash {
    my $matches = {};
    
    $matches->{'transposable_element'}{'dna_transposon'} 
        = [{'Mariner/Tc1' => []}, {'hAT' => []}, 
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
    
    $matches->{'transposable_element'}{'ltr_retrotransposon'} 
        = [{'Gypsy' => []}, {'Copia' => []}, 
	   {'BEL' => []}, {'DIRS' => []}];
    
    $matches->{'transposable_element'}{'endogenous_retrovirus'} 
        = [{'ERV1' => []}, {'ERV2' => []}, 
	   {'ERV3' => []}, {'Lentivirus' => []}, 
	   {'ERV4' => []}];
    
    $matches->{'transposable_element'}{'non-ltr_retrotransposon'} 
        = [{'SINE1/7SL' => []}, {'SINE2/tRNA' => []},
	   {'SINE3/5S' => []},{'SINE4' => []},
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
    
    $matches->{'simple_repeat'}{'Satellite'} 
        = [{'SAT' => []}, {'MSAT' => []}];

    $matches->{'pseudogene'} 
        = [{'rRNA' => []}, {'tRNA' => []}, {'snRNA' => []}];

    $matches->{'integrated_virus'} 
        = [{'DNA_Virus' => []}, {'Caulimoviridae' => []}];
    
    return $matches;
}
