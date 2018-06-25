#!/usr/bin/env perl

# Purpose: Take a GFF3 and FASTA from Tephra (or any source that confines to the Sequence Ontology) 
# and compare them to ensure the number of elements and the IDs in each file are the same. 
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
#dd $fas_ids and exit;

my $re = qr/helitron|(?:non_)?(?:LARD|LTR|TRIM)_retrotransposon|terminal_inverted_repeat_element|MITE/;

my $gff_ids = {};
my $id;
for my $rep_region (keys %$features) {
    my ($chr, $rreg_id, $rreg_start, $rreg_end) = split /\|\|/, $rep_region;
    if ($rreg_id =~ /fragment/) {
	$id = join "_", $rreg_id, $chr, $rreg_start, $rreg_end;
    }
    else {
	for my $feature (@{$features->{$rep_region}}) {
	    if ($feature->{type} =~ /$re/) { 
		my $region = @{$feature->{attributes}{ID}}[0];
		my ($seq_id, $start, $end) = @{$feature}{qw(seq_id start end)};
		if (defined $feature->{attributes}{family}) {
		    my $family = @{$feature->{attributes}{family}}[0];
		    $id = join "_", $family, $region, $seq_id, $start, $end;
		}
		else {
		    $id = join "_", $region, $seq_id, $start, $end;
		}
	    }	
	}
    }
    $gff_ids->{$id} = 1;
}
#dd $gff_ids and exit;

if (%$fas_ids && %$gff_ids) {
    say join "\n", scalar(keys %$fas_ids), scalar(keys %$gff_ids);
    use Test::More tests => 1;
    is_deeply($fas_ids, $gff_ids, 'FASTA and GFF3 IDs are the same');
}

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

    my ($in, $gffio, $header);
    if ($gff =~ /\.gz/) {
        open $in, '-|', 'zcat', $gff or die "\nERROR: Could not open file: $gff\n";
    }
    else {
        open $in, '<', $gff or die "\nERROR: Could not open file: $gff\n";
    }

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
    close $in;
    chomp $header;

    if ($gff =~ /\.gz/) {
        open $gffio, '-|', 'zcat', $gff or die "\nERROR: Could not open file: $gff\n";
    }
    else {
        open $gffio, '<', $gff or die "\nERROR: Could not open file: $gff\n";
    }

    my ($start, $end, $region, $key, %features);
    while (my $line = <$gffio>) {
        chomp $line;
        next if $line =~ /^#/;
	my $feature = gff3_parse_feature( $line );
	#next if $feature->{type} =~ /solo_LTR|similarity/;
	if ($feature->{type} =~ /helitron|non_LTR_retrotransposon|similarity|solo_LTR/) { 
	    ($start, $end) = @{$feature}{qw(start end)};
	    $region = @{$feature->{attributes}{ID}}[0];
	    $key = join "||", $feature->{seq_id}, $region, $start, $end;
	    push @{$features{$key}}, $feature;
	    next;
	}
        if ($feature->{type} eq 'repeat_region') {
            $region = @{$feature->{attributes}{ID}}[0];
            ($start, $end) = @{$feature}{qw(start end)};
	    $key = join "||", $feature->{seq_id}, $region, $start, $end;

        }
	if ($feature->{type} !~ /repeat_region|helitron|non_LTR_retrotransposon|similarity|solo_LTR/) {
	    #dd $feature and exit unless defined $start and defined $end;
            if ($feature->{start} >= $start && $feature->{end} <= $end) {
		push @{$features{$key}}, $feature;
            }
        }
    }
    close $gffio;

    return ($header, \%features);
}
