#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;
use autodie;
use File::Basename;
use Sort::Naturally;
use Bio::GFF3::LowLevel qw(gff3_parse_feature gff3_format_feature);
use Data::Dump::Color;
use Getopt::Long;
use Data::Printer;

my $usage = "\nUSAGE: ".basename($0)." -g gff -l list > filtered.gff3

The output is written to STDOUT, so redirect the output to a file as shown
above.

";

my %opts;
GetOptions(\%opts, 'gff|g=s', 'list|l=s', 'allfeatures|a');

say $usage and exit(1) unless %opts;

open my $in, '<', $opts{list};
my %telist = map { chomp; $_ => 1 } <$in>;
close $in;

my ($header, $features) = collect_gff_features($opts{gff}, \%telist);
#dd $features and exit;
say $header;

for my $chr_id (nsort keys %$features) {
    for my $rep_region (keys %{$features->{$chr_id}}) {
	for my $feature (@{$features->{$chr_id}{$rep_region}}) {

	    my $teid = $feature->{attributes}{ID}[0];
	    next unless defined $teid;

	    if ($teid =~ /(?:LTR|TRIM|LARD)_retrotransposon|terminal_inverted_repeat_element|MITE/) {
		if (exists $telist{$teid}) {
		    _write_features({ parent => $rep_region, features => $features->{$chr_id}{$rep_region} });
		}
	    }
	}
    }
}

sub collect_gff_features {
    my ($gff) = @_;

    my $header;
    open my $in, '<', $gff or die "\n[ERROR]: Could not open file: $gff\n";
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

    open my $gffio, '<', $gff or die "\n[ERROR]: Could not open file: $gff\n";

    my ($seq_id, $start, $end, $phase, $score, $source, $strand, $type, $region, $key, %features);
    while (my $line = <$gffio>) {
        chomp $line;
        next if $line =~ /^#/;
	my $feature = gff3_parse_feature( $line );
	#dd $feature and exit;
        if ($feature->{type} eq 'repeat_region') {
            $region = @{$feature->{attributes}{ID}}[0];
	    #{
		#attributes => { ID => ["repeat_region1"] }, # {0}
		#end        => "7942",                       # {1}
		#phase      => undef,                        # {2}
		#score      => undef,                        # {3}
		#seq_id     => "Ha412HOChr00c00215",         # {4}
		#source     => "LTRharvest",                 # {5}
		#start      => "2653",                       # {6}
		#strand     => "+",                          # {7}
		#type       => "repeat_region",              # {8}
	    #}

            ($seq_id, $start, $end, $phase, $score, $source, $strand, $type) = 
		@{$feature}{qw(seq_id start end phase score source strand type)};
	    $phase //= '.';
	    $score //= '.';
	    
	    $key = join "||", $seq_id, $source, $type, $start, $end, $score, $strand, $phase, $region;
        }
	if ($feature->{type} ne 'repeat_region') {
            if ($feature->{start} >= $start && $feature->{end} <= $end) {
		push @{$features{$seq_id}{$key}}, $feature;
            }
        }
    }
    close $gffio;

    return ($header, \%features);
}

sub _write_features {
    my ($ref) = @_;

    my @rep_features = split /\|\|/, $ref->{parent};
    $rep_features[8] =~ s/^/ID=/;
    my $parent_feat = join "\t", @rep_features;
    chomp $parent_feat;
    say $parent_feat;

    for my $feat (@{$ref->{features}}) {
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
