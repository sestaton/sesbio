#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;
#use autodie;
use Statistics::Descriptive;
use Sort::Naturally;
use File::Spec;
use File::Find;
use File::Basename;
use Bio::DB::HTS::Kseq;
use Bio::DB::HTS::Faidx;
use Bio::GFF3::LowLevel qw(gff3_parse_feature gff3_format_feature);
use List::MoreUtils     qw(indexes any);
use List::Util          qw(min max);
use List::UtilsBy       qw(nsort_by);
use Time::HiRes         qw(gettimeofday);
use File::Path          qw(make_path);
use Cwd;
use Data::Dump::Color;

my $usage = "$0 sf_gff\n";
my $infile = shift or die $usage;

#my $type = 'LTR_retrotransposon';

my ($header, $features) = collect_gff_features($infile);
#dd $features and exit;

my $refct = scalar(keys%$features);
my @remove;
for my $rep_region (keys %$features) {
    for my $feature (@{$features->{$rep_region}}) {
	#if ($feature->{type} =~ /LTR_retrotransposon|solo_LTR/) { # 'terminal_inverted_repeat_element') {
	if ($feature->{type} eq 'helitron') { 
	    push @remove, $rep_region;
	}
    }
}
delete $features->{$_} for @remove;
my $filtct = scalar(keys %$features);
say STDERR ($refct-$filtct),' elements filtered.';

#dd $features and exit;
say $header;
my ($seq_id, $source, $start, $end, $feats, $strand);
for my $rep_region (nsort_by { m/(?:repeat_region)?(?:SO:)?\d+\|\|(\d+)\|\|\d+/ and $1 } keys %$features) {
    if ($rep_region =~ /^SO/) {
	my $feature = shift @{$features->{$rep_region}}; #) {
	my $gff3_str = gff3_format_feature($feature);
	say $gff3_str;
	next;
    }

    my ($rreg_id, $rreg_start, $rreg_end) = split /\|\|/, $rep_region; 
    for my $feature (@{$features->{$rep_region}}) {
	($seq_id, $source, $start, $end, $strand) 
	    = @{$feature}{qw(seq_id source start end strand)};

	my $gff3_str = gff3_format_feature($feature);
	$feats .= $gff3_str;
    }
    chomp $feats;
    say join "\t", $seq_id, $source, 'repeat_region', $rreg_start, $rreg_end, '.', $strand, '.', "ID=$rreg_id";
    say $feats;

    undef $feats;
}

sub collect_gff_features {
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
	    #$region = @{$feature->{attributes}{ID}}[0];
	    $region = @{$feature->{attributes}{Ontology_term}}[0];
	    $key = join "||", $region, $start, $end;
	    push @{$features{$key}}, $feature;
	    next;
	}
        if ($feature->{type} eq 'repeat_region') {
            $region = @{$feature->{attributes}{ID}}[0];
            ($start, $end) = @{$feature}{qw(start end)};
	    $key = join "||", $region, $start, $end;

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

