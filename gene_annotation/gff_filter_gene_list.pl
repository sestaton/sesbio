#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;
use Sort::Naturally;
use Bio::GFF3::LowLevel qw(gff3_parse_feature gff3_format_feature);
#use Data::Dump::Color;

my $usage = "$0 gff list";
my $gff   = shift or die $usage;
my $list  = shift or die $usage;

open my $in, '<', $list;
my %genelist = map { chomp; $_ => 1 } <$in>;
close $in;

my ($header, $features) = collect_gff_features($gff, \%genelist);
#dd $features and exit;
say $header;

for my $chr (nsort keys %$features) {
    for my $id ( map $_->[0],
		 sort { $a->[2] <=> $b->[2] || $a->[1] cmp $b->[1] }
		 map [ $_, split /\|\|/ ],
		 keys %{$features->{$chr}} ) { 

	my $geneid = (split /\|\|/, $id)[0];
	if ($geneid =~ /gene/) {
	    _write_features($features->{$chr}{$id})
		if exists $genelist{$geneid};
	}
	else {
	    _write_features($features->{$chr}{$id});
	}
    }
}

sub collect_gff_features {
    my ($gff, $genelist) = @_;

    my $header;
    open my $in, '<', $gff or die "\nERROR: Could not open file: $gff\n";
    while (<$in>) {
	chomp;
	next if /^###$/;
	if (/^##?\w+/) {
	    $header .= $_."\n";
	}
	else {
	    last;
	}
    }
    close $in;
    chomp $header;

    open my $gffio, '<', $gff or die "\nERROR: Could not open file: $gff\n";

    my $gct = 0;
    my ($geneid, $start, $end, $region, $key, $hash, $type, $source, 
	$seqid, $strand, $score, $attributes, %features);
    while (my $line = <$gffio>) {
        chomp $line;
        next if $line =~ /^#/;
	last if $line =~ /^>/;
	my $feature = gff3_parse_feature( $line );
	if (defined $feature->{type} && 
	    $feature->{type} =~ /gene|expressed_sequence_match|protein_match|^match$/) {
	    $gct++;
	    $geneid = @{$feature->{attributes}{ID}}[0];
	    
	    ($start, $end, $type, $source, $seqid, $strand, $score) =
                @{$feature}{qw(start end type source seq_id strand score)};
            $score //= '.';
            $key = join "||", $geneid, $start, $end, $type, $source, $seqid, $strand, $score;
            $features{$seqid}{$key}{parent} = $feature;
	}
	if (defined $feature->{type} && 
	    $feature->{type} =~ /CDS|exon|five_prime_UTR|match_part|mRNA|three_prime_UTR/) {
	    if (defined $start && $feature->{start} >= $start && 
		defined $end   && $feature->{end}   <= $end) {
		push @{$features{$seqid}{$key}{parts}}, $feature;
	    }
	}
    }
    close $gffio;

    return ($header, \%features);
}

sub _write_features {
    my ($ref) = @_;
    my $parent_feat = gff3_format_feature($ref->{parent});
    chomp $parent_feat;
    say $parent_feat;
    
    for my $feat (@{$ref->{parts}}) {
	my $gff3_str = gff3_format_feature($feat);
	chomp $gff3_str;
	say $gff3_str;
    }
    say '###';
}
