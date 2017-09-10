#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;
use Sort::Naturally;
use Bio::DB::HTS::Kseq;
use Bio::GFF3::LowLevel qw(gff3_parse_feature gff3_format_feature);
use Tephra::Annotation::Util;
#use Data::Dump::Color;

my $usage = "$0 genome sf_gff\n";
my $seqfile = shift or die $usage;
my $gfffile = shift or die $usage;

my ($header, $features) = collect_gff_features($gfffile);

my $util = Tephra::Annotation::Util->new;

my %sfams;
my ($seq_id, $source, $fstart, $fend, $sfam);
for my $chr (nsort keys %$features) {
    for my $rreg (
	map  { $_->[0] }
	sort { $a->[1] <=> $b->[1] }
	map  { [ $_, /\w+(?:\d+)?\|\|(\d+)\|\|\d+/ ] } 
	keys %{$features->{$chr}} ) {
	    
	if ($rreg =~ /helitron|non_LTR_retrotransposon/i) {
	    #dd $features->{$chr} and exit;
	    my $feat = @{$features->{$chr}{$rreg}}[0];
	    #dd $feat and exit;
	    my $fam_id = @{$feat->{attributes}{family}}[0];

	    my $sf;
	    if (defined $fam_id && $fam_id !~ /^0$/) { 
		$sf = $util->map_superfamily_name($fam_id);
	    }
	    elsif ($rreg =~ /helitron/) {
		$sf = 'Helitron';
	    }
	    else {
		$sf = 'LINE';
	    }

	    push @{$sfams{$sf}{$chr}}, $feat;
	}
	elsif ($rreg =~ /repeat_region/) {
	    for my $feature (@{$features->{$chr}{$rreg}}) {
		if ($feature->{type} =~ /^LTR_retrotransposon|TRIM_retrotransposon|terminal_inverted_repeat_element/) { 
		    my $id = defined @{$feature->{attributes}{family}}[0] ? @{$feature->{attributes}{family}}[0] 
		    : @{$feature->{attributes}{ID}}[0];
		    $sfam = $util->map_superfamily_name($id);
		    push @{$sfams{$sfam}{$chr}}, $features->{$chr}{$rreg};
		}
	    }
	}
    }	
}
#dd \%sfams and exit;

for my $sfam (nsort keys %sfams) {
    my $sfam_id = $sfam;
    $sfam_id =~ s/\//-/;
    my $outfile = $sfam_id.'_tephra_features.gff3';
    open my $out, '>', $outfile or die $!;
    say $out $header;

    for my $chr (nsort keys %{$sfams{$sfam}}) {
	for my $elements (@{$sfams{$sfam}{$chr}}) {
	    if (ref($elements) ne 'ARRAY') {
		#dd $elements and exit;
		my $gff_feats = gff3_format_feature($elements);
		chomp $gff_feats;
		say $out $gff_feats;
	    }
	    else {
		for my $feat (@$elements) {
		    my $gff_feats = gff3_format_feature($feat);
		    chomp $gff_feats;
		    say $out $gff_feats;
		}
	    }
	}
    }
    close $out;
}

exit;
#
# methods
#
sub collect_gff_features {
    my ($gff) = @_;

    my $header;
    open my $in, '<', $gff or die "\nERROR: Could not open file: $gff\n";
    #open my $in, '-|', 'zcat', $gff or die $!;

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
    #close $in or $? != 0 or die "close: $!";
    chomp $header;

    open my $gffio, '<', $gff or die "\nERROR: Could not open file: $gff\n";
    #open my $gffio, '-|', 'zcat', $gff or die $!;

    my ($start, $end, $region, $key, %features);
    while (my $line = <$gffio>) {
        chomp $line;
        next if $line =~ /^#/;
	my $feature = gff3_parse_feature( $line );
	next if $feature->{type} =~ /similarity|gene/;
	if ($feature->{type} =~ /helitron|non_LTR_retrotransposon/) { 
	    $region = @{$feature->{attributes}{ID}}[0];
	    $key = join "||", $region, $feature->{start}, $feature->{end};
	    push @{$features{$feature->{seq_id}}{$key}}, $feature;
	    next;
	}
        if ($feature->{type} eq 'repeat_region') {
            $region = @{$feature->{attributes}{ID}}[0];
            ($start, $end) = @{$feature}{qw(start end)};
	    $key = join "||", $region, $start, $end;

        }
	if ($feature->{type} !~ /repeat_region|gene|exon|intron|_utr|cds|rna|similarity|helitron|non_LTR_retrotransposon/i) {
            if ($feature->{start} >= $start && $feature->{end} <= $end) {
		push @{$features{$feature->{seq_id}}{$key}}, $feature;
            }
        }
    }
    close $gffio;

    return ($header, \%features);
}
