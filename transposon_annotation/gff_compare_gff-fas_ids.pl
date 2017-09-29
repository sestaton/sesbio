#!/usr/bin/env perl

# Take a GFF3 and FASTA from Tephra and compare them to ensure the number of elements
# and the IDs in each file are the same. 
#
# TODO: Add checks for length, overlap 

use 5.010;
use strict;
use warnings;
use File::Basename;
use Sort::Naturally;
use Bio::DB::HTS::Kseq;
use Bio::DB::HTS::Faidx;
use Bio::GFF3::LowLevel qw(gff3_parse_feature gff3_format_feature);
#use Data::Dump::Color;
use Carp 'croak';

my $usage = "USAGE: ".basename($0)." gff.gz fasta.gz";
my $gff = shift or die $usage;
my $fas = shift or die $usage;

my ($header, $features) = collect_all_gff_features($gff);
#dd $features and exit;
my $fas_ids = get_ids($fas);

my $gff_ids = {};
for my $rep_region (map  { $_->[0] }
		    sort { ncmp($a->[1], $b->[1]) || $a->[2] <=> $b->[2] }
		    map  { [ $_, /(\w+(?:\d+\w+)?(?:\d+)?)\|\|\w+\d+\|\|(\d+)\|\|\d+/ ] } keys %$features) {

    my ($chr, $rreg_id, $rreg_start, $rreg_end) = split /\|\|/, $rep_region;
    my ($feats, $source, $region, $seq_id, $start, $end, $strand, $type);

    if ($rreg_id =~ /DHH/) {
        my $feature = shift @{$features->{$rep_region}};
	my $region = @{$feature->{attributes}{ID}}[0];
	my ($seq_id, $start, $end, $strand) = @{$feature}{qw(seq_id start end strand)}; 
	my $id = join "_", $region, $seq_id, $start, $end;
	$gff_ids->{$id} = 1;
    }
    else { 

	for my $feature (@{$features->{$rep_region}}) {
	    ($seq_id, $source, $start, $end, $strand)
		= @{$feature}{qw(seq_id source start end strand)};
	    
	    if ($feature->{type} eq 'LTR_retrotransposon') {
		$region = @{$feature->{attributes}{ID}}[0];
		($seq_id, $start, $end) = @{$feature}{qw(seq_id start end)};
		my $family = @{$feature->{attributes}{family}}[0];
		my $id = join "_", $family, $region, $seq_id, $start, $end;
		$gff_ids->{$id} = 1;
	    }
	    
	    if ($feature->{type} eq 'non_LTR_retrotransposon') {
		$region = @{$feature->{attributes}{ID}}[0];
		my $superfamily = @{$feature->{attributes}{Name}}[0];
		($seq_id, $start, $end) = @{$feature}{qw(seq_id start end)};
		my $id = join "_", $superfamily, $region, $seq_id, $start, $end;
		$gff_ids->{$id} = 1;
	    }
	    
	    if ($feature->{type} eq 'terminal_inverted_repeat_element') {
		#HanXRQCP TIRvish terminal_inverted_repeat_element 62679 64963 . ? . ID=terminal_inverted_repeat_element1;Parent=repeat_region1;superfamily=DTX;tir_similarity=89.66
		$region = @{$feature->{attributes}{ID}}[0];
		($seq_id, $start, $end, $strand) = @{$feature}{qw(seq_id start end strand)};
		my $superfamily = @{$feature->{attributes}{superfamily}}[0];
		my $id = join "_", $superfamily, $region, $seq_id, $start, $end;
		$gff_ids->{$id} = 1;
	    }
	}
    }
}

say join "\n", scalar(keys %$fas_ids), scalar(keys %$gff_ids);
use Test::More tests => 1;
is_deeply($fas_ids, $gff_ids, 'FASTA and GFF3 IDs are the same');

exit;
#
# methods
#
sub get_ids {
    my ($fas) = @_;

    my %ids;
    my $kseq = Bio::DB::HTS::Kseq->new($fas);
    my $iter = $kseq->iterator;

    while (my $seqobj = $iter->next_seq) {
	$ids{ $seqobj->{name} } = 1;
    }

    return \%ids;
}

sub collect_all_gff_features {
    my ($gff) = @_;
    my $header;
    #open my $in, '<', $gff or die "\nERROR: Could not open file: $gff\n";
    open my $in, '-|', 'zcat', $gff or die $!;
    while (my $line = <$in>) {
	chomp $line;
	next if $line =~ /^###$/;
	if ($line =~ /^##?\w+/) {
	    $header .= $line."\n";
	}
	else {
	    last;
	}
    }
    #close $in;
    close $in or $? != 0 or die "close: $!";
    chomp $header;

    #open my $gffio, '<', $gff or die "\nERROR: Could not open file: $gff\n";
    open my $gffio, '-|', 'zcat', $gff or die $!;

    my ($start, $end, $region, $key, %features);
    while (my $line = <$gffio>) {
        chomp $line;
        next if $line =~ /^#/;
	my $feature = gff3_parse_feature( $line );
	if ($feature->{type} eq 'helitron') { 
	    $region = @{$feature->{attributes}{ID}}[0];
	    #$region = @{$feature->{attributes}{Ontology_term}}[0];
	    $key = join "||", $feature->{seq_id}, $region, $start, $end;
	    push @{$features{$key}}, $feature;
	    next;
	}
        if ($feature->{type} eq 'repeat_region') {
            $region = @{$feature->{attributes}{ID}}[0];
            ($start, $end) = @{$feature}{qw(start end)};
	    $key = join "||", $feature->{seq_id}, $region, $start, $end;

        }
	if ($feature->{type} ne 'repeat_region') {
            if ($feature->{start} >= $start && $feature->{end} <= $end) {
		push @{$features{$key}}, $feature;
            }
        }
    }
    close $gffio;

    return ($header, \%features);
}
