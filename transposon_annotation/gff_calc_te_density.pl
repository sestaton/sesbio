#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;
use Sort::Naturally;
use Bio::DB::HTS::Kseq;
use Bio::GFF3::LowLevel qw(gff3_parse_feature gff3_format_feature);
use Tephra::Annotation::Util;
use Data::Dump::Color;

my $usage = "$0 genome sf_gff\n";
my $seqfile = shift or die $usage;
my $gfffile = shift or die $usage;

my $windows = make_windows($seqfile);
my ($header, $features) = collect_gff_features($gfffile);
#dd $features and exit;

my %bins;
my $binct = 0;
my ($seq_id, $source, $fstart, $fend);
for my $chr (nsort keys %$windows) {
    for my $window (nsort @{$windows->{$chr}}) {
        my ($wct, $start, $end) = split /\|\|/, $window;
        $start = $start == 0 ? $start+1 : $start;
	my $regions = $features->{$chr};
	for my $rreg (
	    map  { $_->[0] }
	    sort { $a->[1] <=> $b->[1] }
	    map  { [ $_, /\w+(?:\d+)?\|\|(\d+)\|\|\d+/ ] } 
		  keys %$regions ) {

	    my ($rregid, $rstart, $rend) = split /\|\|/, $rreg;
	    if ($rstart >= $start && $rend <= $end) { 
		for my $feature (@{$regions->{$rreg}}) {
		    if ($feature->{type} =~ /^LTR_retrotransposon/) { 
                        ($seq_id, $source, $fstart, $fend)
                            = @{$feature}{qw(seq_id source start end)};
                        my $len = $fend-$fstart+1;
			push @{$bins{$chr}{$window}{ $feature->{type} }},
			    { length => $len, ID => @{$feature->{attributes}{family}}[0] };
                    }
		    if ($feature->{type} =~ /terminal_inverted_repeat_element/) {
			 ($seq_id, $source, $fstart, $fend)
			     = @{$feature}{qw(seq_id source start end)};
			 my $len = $fend-$fstart+1;
			 push @{$bins{$chr}{$window}{ $feature->{type} }},
			     { length => $len, ID => @{$feature->{attributes}{superfamily}}[0] };
		    }
		    if ($feature->{type} =~ /non_LTR_retrotransposon|TRIM_retrotransposon|helitron|similarity/) {
			($seq_id, $source, $fstart, $fend) 
			    = @{$feature}{qw(seq_id source start end)};
			my $len = $fend-$fstart+1;
			push @{$bins{$chr}{$window}{ $feature->{type} }}, 
			    { length => $len, ID => @{$feature->{attributes}{ID}}[0] };
		    }
		    if ($feature->{type} eq 'gene') {
                        ($seq_id, $source, $fstart, $fend)
                            = @{$feature}{qw(seq_id source start end)};
                        my $len = $fend-$fstart+1;
			if ($len >= 500 && $len <= 20000) {
			    push @{$bins{$chr}{$window}{ $feature->{type} }},
			    { length => $len, ID => @{$feature->{attributes}{ID}}[0] };
			}
		    }
		}
	    }
	}		
    }
}
#dd \%bins and exit;

my $util = Tephra::Annotation::Util->new;
my %reduced;

for my $chr (nsort keys %bins) {
    for my $bin (nsort keys %{$bins{$chr}}) {
	for my $type (keys %{$bins{$chr}{$bin}}) {
	    for my $match (@{$bins{$chr}{$bin}{$type}}) {
		my $id;
		if ($match->{ID} =~ /trim/i) {
		    $id = 'TRIM';
		}
		elsif ($match->{ID} =~ /^gene/) {
		    $id = 'gene';
		}
		else {
		    $id = $util->map_superfamily_name($match->{ID}); 
		}
		$reduced{$id}{$chr}{$bin} += $match->{length};
	    }
	}
    }
}

#dd \%reduced and exit;
#say STDERR join "\n", nsort keys %reduced;

for my $type (nsort keys %reduced) {
    my $file = $type.'_1Mbp_bins_circos_density.tsv';
    $file =~ s/\//-/g;
    open my $out, '>', $file or die "\nERROR: Could not open file: $file\n";;
    for my $chr (nsort keys %{$reduced{$type}}) {
	for my $bin (map  { $_->[0] }
		     sort { $a->[1] <=> $b->[1] }
		     map  { [ $_, /(\d+)\|\|\d+\|\|\d+/ ] } keys %{$reduced{$type}{$chr}}) {
	    my ($bin_num, $bin_s, $bin_e) = split /\|\|/, $bin;
	    $bin_s = $bin_s > 0 ? $bin_s : 1;
	    my $bin_length = $bin_e - $bin_s + 1;
	    my $perc = sprintf("%.6f", $reduced{$type}{$chr}{$bin}/$bin_length);
	    say $out join "\t", $chr, $bin_s, $bin_e, $perc;
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
	if ($feature->{type} =~ /helitron|similarity|gene/) { 
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
	if ($feature->{type} !~ /repeat_region|gene|exon|intron|_utr|cds|rna|similarity|helitron/i) {
            if ($feature->{start} >= $start && $feature->{end} <= $end) {
		push @{$features{$feature->{seq_id}}{$key}}, $feature;
            }
        }
    }
    close $gffio;

    return ($header, \%features);
}

sub make_windows {
    my ($sequence) = @_;

    my $wsize = 1e6;

    my %windows;
    my $kseq = Bio::DB::HTS::Kseq->new($sequence);
    my $iter = $kseq->iterator;

    while (my $seqobj = $iter->next_seq) {
	my $seq = $seqobj->{seq};
	my $id  = $seqobj->{name};
	my $length = length($seq);
	my $remainder = $length;

        my ($total, $start, $end, $wct, $chunk_size) = (0, 0, 0, 0, 0);
        my $steps = sprintf("%.0f", $length/$wsize);
	$steps = $length % $wsize == 0 ? $steps : $steps + 1;
        $steps = 0 if $length <= $wsize;

	for (0..$steps) {
	    last if $remainder == 0;
            if ($remainder < $wsize) {
		$start = $end;
		$end = $length;
		$remainder = 0;
	    }
	    else {
		if ($wct) {
		    $start = $end;
		    $end = $start + $wsize;
		    $remainder -= $wsize;
		}
		else {
		    $start = 0;
		    $end = $wsize;
		    $remainder -= $wsize;
		}
	    }
	    push @{$windows{ $seqobj->{name} }}, join "||", $wct, $start, $end;
	    $wct++;
	}
    }
    
    return \%windows;
}
