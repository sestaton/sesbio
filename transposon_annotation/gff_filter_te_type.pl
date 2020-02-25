#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;
#use autodie;
use Sort::Naturally;
use Bio::GFF3::LowLevel qw(gff3_parse_feature gff3_format_feature);
use List::UtilsBy       qw(nsort_by);
#use Data::Dump::Color;

my $usage = "$0 sf_gff te_type\n";
my $infile = shift or die $usage;
my $type = shift or die $usage; ## TODO: Check against hash of known SO terms

my ($header, $features) = collect_gff_features($infile);
#dd $features and exit;

my ($refct, $filtct) = (0, 0);

my @remove;
for my $chr (keys %$features) { 
    for my $rep_region (keys %{$features->{$chr}}) {
	$refct++;
	for my $feature (@{$features->{$chr}{$rep_region}}) {
	    if ($feature->{type} eq $type) { 
		push @remove, join "~~", $chr, $rep_region;
		$filtct++;
	    }
	}
    }
}

for my $to_remove (@remove) {
    my ($chr, $rep_region) = split /\~\~/, $to_remove;
    delete $features->{$chr}{$rep_region};
}

my $remct = @remove;
#my $filtct = scalar(keys %$features);
say STDERR "$refct reference features; $filtct elements identified for removal";
say STDERR ($refct-$filtct),' elements filtered.';
say STDERR "$remct removed";

#dd $features and exit;
say $header;

my ($id, $seq_id, $source, $start, $end, $feats, $strand);
for my $chr (nsort keys %$features) { 
    for my $rep_region (
	map  { $_->[0] }
	sort { $a->[1] <=> $b->[1] }
	map  { [ $_, /\S+\|\|(\d+)\|\|\d+/ ] } 
	keys %{$features->{$chr}}) {

	my ($rreg_id, $rreg_start, $rreg_end) = split /\|\|/, $rep_region;
	if ($rreg_id =~ /helitron|non_LTR_retrotransposon|fragment|solo_LTR/) {
	    $id = join "_", $rreg_id, $chr, $rreg_start, $rreg_end;
	    my $feature = shift @{$features->{$chr}{$rep_region}};
	    my $gff3_str = gff3_format_feature($feature);
	    chomp $gff3_str;
	    say $gff3_str;
	    say '##';
	}
	else {
    	    for my $feature (@{$features->{$chr}{$rep_region}}) {
		($seq_id, $source, $start, $end, $strand) 
		    = @{$feature}{qw(seq_id source start end strand)};
		
		my $gff3_str = gff3_format_feature($feature);
		$feats .= $gff3_str;
	    }
	    chomp $feats;
	    say join "\t", $seq_id, $source, 'repeat_region', $rreg_start, $rreg_end, '.', $strand, '.', "ID=$rreg_id";
	    say $feats;
	    say '##';
	    undef $feats;
	}
    }
}

exit;

#
# methods
#
sub collect_gff_features {
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
    #close $in;
    close $in or $? != 0 or die "close: $!";
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
	if ($feature->{type} =~ /helitron|non_LTR_retrotransposon|similarity|solo_LTR/) { 
	    ($start, $end) = @{$feature}{qw(start end)};
	    $region = @{$feature->{attributes}{ID}}[0];
	    $key = join "||", $region, $start, $end;
	    push @{$features{ $feature->{seq_id} }{ $key }}, $feature;
	    next;
	}
        if ($feature->{type} eq 'repeat_region') {
            $region = @{$feature->{attributes}{ID}}[0];
            ($start, $end) = @{$feature}{qw(start end)};
	    $key = join "||", $region, $start, $end;

        }
	#if ($feature->{type} ne 'repeat_region') {
	if ($feature->{type} !~ /repeat_region|helitron|non_LTR_retrotransposon|similarity|solo_LTR/) {
            if ($feature->{start} >= $start && $feature->{end} <= $end) {
		push @{$features{ $feature->{seq_id} }{ $key }}, $feature;
            }
        }
    }
    close $gffio;

    return ($header, \%features);
}

