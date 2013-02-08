#!/usr/bin/env perl

use 5.014;
use strict;
use warnings;
use autodie qw(open);
use lib qw(/home/jmblab/statonse/apps/perlmod/Data-Dump-1.21/blib/lib);
use Data::Dump qw(dd);
use JSON;

my $usage = "$0 idlist\n";
my $infile = shift or die $usage;
open(my $in, '<', $infile);

my $matches = build_repbase_hash();

#dd $matches; # correct here
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
close($in);

#dd %family_map;

for my $type (keys %$matches) {
    #say $type;  ## correct here
    unless ($type eq 'pseudogene' || $type eq 'integrated_virus') {
	for my $class (keys %{$matches->{$type}}) {
	    #say $class;
	    for my $superfam (@{$matches->{$type}{$class}}) {
		#say $superfam;
		for my $superfam_h (keys %$superfam) {
		    #say $superfam_h;
		    my $superfam_cp = lc($superfam_h);

		    # begin walking through map 
		    #
		    for my $mapped_fam (keys %family_map) {
			
			my $mapped_fam_cp = lc($mapped_fam);
			#$mapped_fam_cp =~ s/\d\///;
			#say 'Here is superfam_cp: ',$superfam_cp,' Here is mapped_fam_cp: ',$mapped_fam_cp;

			if (length($superfam_cp) > 1 && length($mapped_fam_cp) > 1) {
			    
			    if ($mapped_fam_cp =~ /sine/i && $superfam_cp =~ /sine/i) {
				say 'Here is superfam_cp: ',$superfam_cp,' Here is mapped_fam_cp: ',$mapped_fam_cp; 
				
				#for my $sine_fam (keys %$superfam_h) {
				#	for my $sine_fam_mem (@{$superfam_h->{$sine_fam}}) {
				
				#	    if (lc($sine_fam_mem) =~ lc($mapped_fam)) {
				#push @{$matches->{$type}{$class}{$superfam}{$sine_fam}{$sine_fam_mem}}, $family_map{$mapped_fam};
				#push @{$matches->{$type}{$sine_fam}}, $family_map{$mapped_fam};
				#		push @{$superfam_h->{$sine_fam}{$sine_fam_mem}}, $family_map{$mapped_fam};
				#	    }
				#	    else {
				#		say "WARNING: ",lc($sine_fam_mem)," is not matching ",lc($mapped_fam);
				#	    }
				#	} 
				#    } 
			    }
			    #else {
			    elsif ($mapped_fam_cp =~ /$superfam_cp/) {
				say 'Here is superfam_cp: ',$superfam_cp,' Here is mapped_cp: ',$mapped_fam_cp;
				#for my $nonsine_fam (keys %$superfam_h) {
				#	if (lc($nonsine_fam) =~ lc($mapped_fam)) {
				#push @{$matches->{$type}{$class}{$superfam}{$sf_h}}, $family_map{$mapped_fam};
				#push @{$superfam_fam->{$fam_h}}, $family_map{$mapped_fam};
				#	    push @{$superfam_h->{$nonsine_fam}}, $family_map{$mapped_fam};
				#	}
				#	else {
				#	    say "WARNING2: ",lc($nonsine_fam)," is not matching ",lc($mapped_fam);
				#	}
			    }
			}
			elsif (length($mapped_fam_cp) == 1 && length($superfam_cp) == 1) {
			    if ($mapped_fam_cp =~ /$superfam_cp/) {
                                say 'Here is superfam_cp: ',$superfam_cp,' Here is mapped_cp: ',$mapped_fam_cp;
			    }
			}
		    }
		}
	    }
	} 
    } 
}
#dd $matches;
    
#my $json = JSON->new->utf8->space_after->encode($matches);
#my $json = encode_json $matches;
#say $json;

#my $hash = JSON->new->utf8->space_after->decode($json);
#dd $hash;

#
# subs
#
sub build_repbase_hash {

    my $matches = {};
    
    $matches->{'transposable_element'}{'dna_transposon'} = [{'Mariner/Tc1' => []}, {'hAT' => []}, {'MuDR' => []}, {'EnSpm' => []}, {'piggyBac' => []}, {'P' => []}, {'Merlin' => []}, {'Harbinger' => []}, {'Transib' => []}, {'Novosib' => []}, {'Helitron' => []}, {'Polinton' => []}, {'Kolobok' => []}, {'ISL2EU' => []}, {'Crypton' => []}, {'Sola' => []}, {'Zator' => []}, {'Ginger/1' => []}, {'Ginger2/TDD' => []}, {'Academ' => []}, {'Zisupton' => []}, {'IS3EU' => []}];
    
    $matches->{'transposable_element'}{'ltr_retrotransposon'} = [{'Gypsy' => []}, {'Copia' => []}, {'BEL' => []}, {'DIRS' => []}];
    
    $matches->{'transposable_element'}{'endogenous_retrovirus'} = [{'ERV1' => []}, {'ERV2' => []}, {'ERV3' => []}, {'Lentivirus' => []}, {'ERV4' => []}];
    
    $matches->{'transposable_element'}{'non-ltr_retrotransposon'} = [{'SINE' => [{'SINE1/7SL' => []},{'SINE2/tRNA' => []},{'SINE3/5S' => []},{'SINE4' => []}]},{'CRE' => []}, {'NeSL' => []}, {'R4' => []}, {'R2' => []}, {'L1' => []}, {'RTE' => []}, {'I' => []}, {'Jockey' => []}, {'CR1' => []}, {'Rex1' => []}, {'RandI' => []}, {'Penelope' => []}, {'Tx1' => []}, {'RTEX' => []}, {'Crack' => []}, {'Nimb' => []}, {'Proto1' => []}, {'Proto2' => []}, {'RTETP' => []}, {'Hero' => []}, {'L2' => []}, {'Tad1' => []}, {'Loa' => []}, {'Ingi' => []}, {'Outcast' => []}, {'R1' => []}, {'Daphne' => []}, {'L2A' => []}, {'L2B' => []}, {'Ambal' => []}, {'Vingi' => []}, {'Kiri' => []}];

    $matches->{'simple_repeat'}{'Satellite'} = [{'SAT' => []}, {'MSAT' => []}];

    $matches->{'pseudogene'} = [{'rRNA' => []}, {'tRNA' => []}, {'snRNA' => []}];

    $matches->{'integrated_virus'} = [{'DNA_Virus' => []}, {'Caulimoviridae' => []}];
    
    return $matches;
}
