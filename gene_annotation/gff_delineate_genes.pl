#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;
use autodie;
use File::Basename;
use Sort::Naturally;
use Bio::GFF3::LowLevel qw(gff3_parse_feature gff3_format_feature);
#use Data::Dump::Color;
use Getopt::Long;

my $usage = "\nUSAGE: ".basename($0)." -g gff\n";

my %opts;
GetOptions(\%opts, 'gff|g=s'); 

say $usage and exit(1) unless %opts;

my ($header, $features) = collect_gff_features($opts{gff});
#dd $features and exit;
say $header;
#exit;

my $gene;
my ($exon_ct, $cds_ct) = (0, 0);
for my $chr (nsort keys %$features) {
    for my $id ( map $_->[0],
		 sort { $a->[2] <=> $b->[2] || $a->[1] cmp $b->[1] }
		 map [ $_, split /\|\|/ ],
		 keys %{$features->{$chr}} ) { 
	
	my $parent_feat = gff3_format_feature($features->{$chr}{$id}{parent});
	chomp $parent_feat;
	say $parent_feat;
    
	for my $feat (@{$features->{$chr}{$id}->{parts}}) {
	    my $gff3_str = gff3_format_feature($feat);
	    chomp $gff3_str;
	    say $gff3_str;
	}
	say '###';

    }
}

exit;
#
# methods
#
sub collect_gff_features {
    my ($gff) = @_;

    my $fh = _get_fh($gff);

    my $header;
    while (<$fh>) {
	chomp;
	next if /^###$/;
	if (/^##?\w+/) {
	    $header .= $_."\n";
	}
	else {
	    last;
	}
    }
    close $fh;
    chomp $header;

    my $gffio = _get_fh($gff);

    my $gct = 0;
    my ($geneid, $start, $end, $region, $key, $hash, $type, $source, 
	$seqid, $strand, $score, $attributes, %features, @remove);

    #my ($exon_ct, $cds_ct) = (0, 0);
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
	    $feature->{type} =~ /CDS|exon|intron|UTR|match|mRNA/) { # matches: five_prime_UTR and three_prime_UTR; 
                                                                    # match and match_part
	    my $attrid = @{$feature->{attributes}{ID}}[0];
	    if (defined $attrid && $attrid =~ /ncrna/i) { # to remove
		push @remove, $key;
	    }

	    #$exon_ct++ if $feature->{type} eq 'exon';
	    #$cds_ct++ if $feature->{type} eq 'CDS';

	    if (defined $start && $feature->{start} >= $start && 
		defined $end   && $feature->{end}   <= $end) {
		push @{$features{$seqid}{$key}{parts}}, $feature;
	    }
	}
	#if (defined $feature->{type} && 
	    #$feature->{type} =~ /ncrna/i) { # to remove
	    #push @remove, $key;
        #}
	
    }
    close $gffio;

    if (@remove) {
	for my $seqid (keys %features) {
	    delete $features{$seqid}{$_} for @remove;
	}
    }

    return ($header, \%features);
}

sub _get_fh {
    my ($gff) = @_;

    my $fh;
    if ($gff =~ /\.gz$/) {
        open $fh, '-|', 'zcat', $gff;
    }
    else {
        open $fh, '<', $gff;
    }

    return $fh;
}
