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

my $usage = "\nUSAGE: ".basename($0)." -g genes.gff3.gz -l list <--allfeatures> <--invert> > filtered.gff3

By default, only gene features are written to the output. With the 
--allfeaures option, all other feautures in the GFF, like BLAST alignments
or other predictions, will be output. The --invert option will write all
genes/features not in <list>.

The output is written to STDOUT, so redirect the output to a file as shown
above.

";

my %opts;
GetOptions(\%opts, 'gff|g=s', 'list|l=s', 'invert|v', 'allfeatures|a');

say $usage and exit(1) unless %opts;

open my $in, '<', $opts{list};
my %genelist = map { chomp; $_ => 1 } <$in>;
close $in;
#dd \%genelist and exit;

my ($header, $features) = collect_gff_features($opts{gff}, \%genelist);
#dd $features and exit;
say $header;

for my $chr (nsort keys %$features) {
    for my $id ( map $_->[0],
		 sort { $a->[2] <=> $b->[2] || $a->[1] cmp $b->[1] }
		 map [ $_, split /\|\|/ ],
		 keys %{$features->{$chr}} ) { 

	my $gene_feat = (split /\|\|/, $id)[2];
	my $geneid = (split /\|\|/, $id)[0];
	if (exists $genelist{$geneid} && !exists $opts{invert}) { 
	    _process_feature_args($features, $chr, $id, $gene_feat, $opts{allfeatures});
	}
	elsif (!exists $genelist{$geneid} && exists $opts{invert}) {
	    _process_feature_args($features, $chr, $id, $gene_feat, $opts{allfeatures});
	}
    }
}

exit;
#
# methods
#
sub collect_gff_features {
    my ($gff, $genelist) = @_;

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

sub _process_feature_args {
    my ($features, $chr, $id, $gene_feat, $allfeatures) = @_; 

    if ($gene_feat eq 'gene') {
	_write_features($features->{$chr}{$id});
    }
    else {
	_write_features($features->{$chr}{$id})
	    if $allfeatures;
    }

    return;
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
