#!/usr/bin/env perl

## NB: This is for testing new annotation features in Transposome,
##     so it is probably not of interest to anyone but me. This is also
##     not likely to be up to date with Transposome, sorry.

use 5.010;
use strict;
use warnings;
use autodie;
use Data::Dump;
use JSON;
use List::MoreUtils qw(first_index);

my $usage   = "$0 fastadb tophit\n";
my $infile  = shift or die $usage;
my @top_hits = qw(BEL1_I_AG BEL2_LTR_AG GYPSY1_LTR_AG Copia_7_AG_LTR MTANGA_I GYPSY32_LTR_AG PegasusA DNA_2_AG Clu_15B_AG);

my ($repeats, $type_map) = _map_repeats($infile);

my (@tops, @annos);
for my $top_hit (@top_hits) {
    my %anno_data = ( filebase     => 'CL100',
		      top_hit      => $top_hit,
		      top_hit_frac => 0.45,
		      readct       => 100000,
		      repeat_map   => $repeats,
		      repeat_type  => $type_map->{$top_hit} );

    my ($top_hit_superfam, $cluster_annot) = blast_to_annotation(\%anno_data);
    push @annos, $cluster_annot;
    push @tops, $top_hit_superfam;
}

dd \@annos;
dd \@tops;

sub blast_to_annotation { 
    my ($anno_data) = @_;
    my $repeats = $anno_data->{repeat_map};

    my $top_hit_superfam = {};
    my $cluster_annot = {};

    keys %$repeats;
    
    for my $type (keys %$repeats) {
	if (defined $anno_data->{repeat_type} && $type eq $anno_data->{repeat_type}) { 
	    if ($type =~ /(pseudogene|integrated_virus)/) {
		$anno_data->{'class'} = $1;
		($top_hit_superfam, $cluster_annot) = _map_hit_family($repeats->{$type}, $anno_data);
	    }
	    elsif ($type =~ /(satellite)/i) {
		$anno_data->{'class'} = $1;
		($top_hit_superfam, $cluster_annot) = _map_hit_family($repeats->{$type}, $anno_data);
	    }
	    elsif ($type eq 'ltr_retrotransposon') {
		$anno_data->{'class'} = $type;
		($top_hit_superfam, $cluster_annot) = _map_hit_family($repeats->{$type}, $anno_data);
	    }
	    elsif ($type eq 'non-ltr_retrotransposon') {
		$anno_data->{'class'} = $type;
		($top_hit_superfam, $cluster_annot) = _map_hit_family($repeats->{$type}, $anno_data);
	    }
	    elsif ($type eq 'endogenous_retrovirus') {
		$anno_data->{'class'} = $type;
		($top_hit_superfam, $cluster_annot) = _map_hit_family($repeats->{$type}, $anno_data);
	    }
	    elsif ($type eq 'dna_transposon') {
		$anno_data->{'class'} = $type;
		($top_hit_superfam, $cluster_annot) = _map_hit_family($repeats->{$type}, $anno_data);
	    }
	    else {
		my $unk_fam = q{ };
		$anno_data->{'class'} = 'unknown';
		($top_hit_superfam, $cluster_annot) = _map_hit_family($repeats->{$type}, $anno_data);
	    }
	}
    }

    return ($top_hit_superfam, $cluster_annot);
}

sub _map_te_type {
    my ($type) = @_;

    my %map = (
	       'ltr_retrotransposon'     => 'transposable_element',
	       'dna_transposon'          => 'transposable_element',
	       'non-ltr_retrotransposon' => 'transposable_element',
	       'endogenous_retrovirus'   => 'transposable_element',
	       'Satellite'               => 'simple_repeat',
	       );

    return $map{$type} if exists $map{$type};
    return $type if !exists $map{$type};
}

sub _map_hit_family {
    my ($arr_ref, $anno_data) = @_;
    my (%top_hit_superfam, %cluster_annot);
    my ($filebase, $top_hit, $top_hit_frac, $readct, $class) 
        = @{$anno_data}{qw(filebase top_hit top_hit_frac readct class)};

    my $type = _map_te_type($class);
    #say STDERR join q{ }, "DEBUG: ", $top_hit, $class, $type;

    if ($class =~ /pseudogene|integrated_virus/) {
        my $anno_val = mk_key($filebase,
			      $type,
			      $class,
			      '-',
			      '-',
			      $top_hit,
			      $top_hit_frac);
	my $anno_key = mk_key($top_hit, $readct);
	$cluster_annot{$anno_key} = $anno_val;
        return (undef, \%cluster_annot);
    }
    else {
	for my $superfam_h (@$arr_ref) {
            for my $superfam (keys %$superfam_h) {
                for my $family (@{$superfam_h->{$superfam}}) {
                    if ($top_hit =~ /$family/) {
			#say join q{ }, $filebase, $top_hit, $type, $class, $superfam;
			my $family_name = _map_family_name($family);
                        $top_hit_superfam{$top_hit} = $superfam;
                        my $anno_key = mk_key($filebase, 
					      $type,
					      $class,
					      $superfam, 
					      $family_name, 
					      $top_hit, 
					      $top_hit_frac);
        
			#say $anno_key;
			#my $anno_key = mk_key($top_hit, $readct);
			$cluster_annot{$readct} = $anno_key;
                    }
                }
            }
        }
        
        return (\%top_hit_superfam, \%cluster_annot);
    }
}

sub mk_key {
    return join "~~", map { $_ // " " } @_;
}

sub mk_vec {
    my $key = @_;
    return split /\~\~/, $key;
}

sub _map_family_name {
    my ($family) = @_;

    my $family_name;
    if ($family =~ /(^RL[GCX][_-][a-zA-Z]*\d*?[_-]?[a-zA-Z-]+?\d*?)/) {
	$family_name = $1;
    }
    elsif ($family =~ /(^D[HT][ACHMT][_-][a-zA-Z]*\d*?)/) {
	$family_name = $1;
    }
    elsif ($family =~ /(^PPP[_-][a-zA-Z]*\d*?)/) {
	$family_name = $1;
    }
    elsif ($family =~ /(^R[IS][LT][_-][a-zA-Z]*\d*?)/) {
	$family_name = $1;
    }
    else {
	$family_name = $family;
    }

    $family_name =~ s/_I// if $family_name =~ /_I_|_I$/;
    $family_name =~ s/_LTR// if $family_name =~ /_LTR_|_LTR$/;

    return $family_name;
}

sub _map_repeats {
    my ($infile) = @_;

    my $matches = build_repbase_hash();
    my $repeats = map_repeat_types($matches);
    #dd $repeats and exit;

    open my $in, '<', $infile;
    my (%family_map, %type_map, %seen);

    while (my $line = <$in>) {
	chomp $line;
	if ($line =~ /^>/) {
	    $line =~ s/>//;
	    my ($f, $sf, $source)  = split /\t/, $line;
	    next unless defined $sf && defined $f;
	    if ($sf =~ /(\s+)/) {
		$sf =~ s/$1/\_/;
	    }
	    $f =~ s/\s/\_/;
	    push @{$family_map{$sf}}, $f;
	}
    }
    close $in;
    
    for my $mapped_sfam (keys %family_map) {
	my $mapped_sfam_cp = lc($mapped_sfam);
	for my $mapped_fam (@{$family_map{$mapped_sfam}}) {
	    for my $class (keys %$repeats) {
		for my $sfamh (@{$repeats->{$class}}) {
		    my $sfam_index = first_index { $_ eq $sfamh } @{$repeats->{$class}};
		    for my $sfamname (keys %$sfamh) {
			if (lc($sfamname) eq $mapped_sfam_cp) {
			    push @{$repeats->{$class}[$sfam_index]{$sfamname}}, $mapped_fam;
			    $type_map{$mapped_fam} = $class;
			    #say join q{ }, "DEBUG: ", $mapped_fam, $class;
			}
			elsif ($class eq $mapped_sfam_cp) {
			    my $unk_idx = first_index { $_ eq 'unclassified' } @{$repeats->{$class}};
			    push @{$repeats->{$class}[$unk_idx]{'unclassified'}}, $mapped_fam
				unless exists $seen{$mapped_fam};
			    $seen{$mapped_fam} = 1;
			    $type_map{$mapped_fam} = $class;
			    #say join q{ }, "Debug1: ", $mapped_fam, $class;
			}
		    }
		}
	    }
	}
    }
    return ($repeats, \%type_map);
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
		    $repeats{'non-ltr_retrotransposon'} = $matches->{$type}{$tes};
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
            $repeats{'pseudogene'} = $matches->{$type};
        }
        elsif ($type eq 'integrated_virus') { 
            $repeats{'integrated_virus'} = $matches->{$type};
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
	   {'Zisupton' => []}, {'IS3EU' => []}, 
	   {'unclassified' => []}];
    
    $matches->{'transposable_element'}{'ltr_retrotransposon'} 
        = [{'Gypsy' => []}, {'Copia' => []}, 
	   {'BEL' => []}, {'DIRS' => []}, {'unclassified' => []}];
    
    $matches->{'transposable_element'}{'endogenous_retrovirus'} 
        = [{'ERV1' => []}, {'ERV2' => []}, 
	   {'ERV3' => []}, {'Lentivirus' => []}, 
	   {'ERV4' => []}, {'unclassified' => []}];
    
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
	   {'Vingi' => []}, {'Kiri' => []}, {'unclassified' => []}];
    
    $matches->{'simple_repeat'}{'Satellite'} 
        = [{'SAT' => []}, {'MSAT' => []}];

    $matches->{'pseudogene'} 
        = [{'rRNA' => []}, {'tRNA' => []}, {'snRNA' => []}];

    $matches->{'integrated_virus'} 
        = [{'DNA_Virus' => []}, {'Caulimoviridae' => []}];
    
    return $matches;
}
