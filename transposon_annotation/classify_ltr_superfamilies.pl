#!/usr/bin/env perl

use 5.020;
use warnings;
use autodie;
use File::Basename;
use Statistics::Descriptive;
use Bio::Tools::GFF;
use List::UtilsBy qw(nsort_by);
use Data::Dump;
use experimental 'signatures';

my $usage = "$0 gff";
my $gff   = shift or die $usage;
my $header;

open my $in, '<', $gff;
while (<$in>) {
    chomp;
    if (/^#/) {
	$header .= $_."\n";
    }
    else {
	last;
    }
}
close $in;
chomp $header;

my $gffio = Bio::Tools::GFF->new( -file => $gff, -gff_version => 3 );

my ($start, $end, $region, %feature);
while (my $feature = $gffio->next_feature()) {
    if ($feature->primary_tag eq 'repeat_region') {
	my @string = split /\t/, $feature->gff_string;
	($region) = ($string[8] =~ /ID=?\s+?(repeat_region\d+)/);
	($start, $end) = ($feature->start, $feature->end);
    }
    next $feature unless defined $start && defined $end;
    if ($feature->primary_tag ne 'repeat_region') {
	if ($feature->start >= $start && $feature->end <= $end) {
	    push @{$feature{$region.".".$start.".".$end}}, join "||", split /\t/, $feature->gff_string;
	}
    }
}

my $all_ct = (keys %feature);
find_gypsy(\%feature, $header, $gff);
my $gyp_ct = (keys %feature);
find_copia(\%feature, $header, $gff);
my $cop_ct = (keys %feature);
write_unclassified_ltrs(\%feature, $header, $gff);
my $rem_ct = (keys %feature);

say STDERR join "\t", "all", "after_gypsy", "after_copia", "after_rem";
say STDERR join "\t", $all_ct, $gyp_ct, $cop_ct, $rem_ct;
#
# methods
#
sub find_gypsy ($feature, $header, $gff) {
    my @lengths;
    my $gyp_feats;
    my $is_gypsy = 0;
    my $has_pdoms  = 0;
    my $pdoms      = 0;
    
    my ($name, $path, $suffix) = fileparse($gff, qr/\.[^.]*/);
    my $outfile = $name."_gypsy.gff3";
    open my $out, '>>', $outfile;
    say $out $header;

    for my $ltr (nsort_by { m/repeat_region(\d+)/ and $1 } keys %$feature) {
	my ($rreg, $s, $e) = split /\./, $ltr;
	my $len = ($e - $s) + 1;
	my $region = @{$feature->{$ltr}}[0];
	my ($loc, $source) = (split /\|\|/, $region)[0..1];

	for my $feat (@{$feature->{$ltr}}) {
	    my @feats = split /\|\|/, $feat;
	    $feats[8] =~ s/\s\;\s/\;/g;
	    $feats[8] =~ s/\s+/=/g;
	    $feats[8] =~ s/\s+$//;
	    $feats[8] =~ s/=$//;
	    $feats[8] =~ s/=\;/;/g;
	    $feats[8] =~ s/\"//g;
	    if ($feats[2] =~ /protein_match/ && $feats[8] =~ /name=RVT_1|name=Chromo/i) {
		$is_gypsy = 1;
		$has_pdoms = 1;
	    }
	    $gyp_feats .= join "\t", @feats, "\n";
	}
	if ($is_gypsy) {
	    chomp $gyp_feats;
	    say $out join "\t", $loc, $source, 'repeat_region', $s, $e, '.', '?', '.', "ID=$rreg";
	    say $out $gyp_feats;
	    delete $feature->{$ltr};
	    push @lengths, $len;
	    $pdoms++ if $has_pdoms;
	}
	undef $gyp_feats;
	$is_gypsy = 0;
	$has_pdoms  = 0;
    }
    close $out;

    say STDERR join "\t", "gypsy_count", "min_length", "max_length", "mean_length", "elements_with_protein_matches";
    my $stat = Statistics::Descriptive::Full->new;
    $stat->add_data(@lengths);
    my $min   = $stat->min;
    my $max   = $stat->max;
    my $mean  = $stat->mean;
    my $count = $stat->count;
    say STDERR join "\t", $count, $min, $max, sprintf("%.2f", $mean), $pdoms;
}

sub find_copia ($feature, $header, $gff) {
    my @lengths;
    my $cop_feats;
    my $is_copia = 0;
    my $has_pdoms  = 0;
    my $pdoms      = 0;
    
    my ($name, $path, $suffix) = fileparse($gff, qr/\.[^.]*/);
    my $outfile = $name."_copia.gff3";
    open my $out, '>>', $outfile;
    say $out $header;

    for my $ltr (nsort_by { m/repeat_region(\d+)/ and $1 } keys %$feature) {
        my ($rreg, $s, $e) = split /\./, $ltr;
        my $len = ($e - $s) + 1;
        my $region = @{$feature->{$ltr}}[0];
        my ($loc, $source) = (split /\|\|/, $region)[0..1];

        for my $feat (@{$feature->{$ltr}}) {
	    my @feats = split /\|\|/, $feat;
	    $feats[8] =~ s/\s\;\s/\;/g;
	    $feats[8] =~ s/\s+/=/g;
	    $feats[8] =~ s/\s+$//;
	    $feats[8] =~ s/=$//;
	    $feats[8] =~ s/=\;/;/g;
	    $feats[8] =~ s/\"//g;
	    if ($feats[2] =~ /protein_match/ && $feats[8] =~ /name=RVT_2/) {
		$is_copia = 1;
		$has_pdoms = 1;
	    }
	    $cop_feats .= join "\t", @feats, "\n";
	}
        if ($is_copia) {
	    chomp $cop_feats;
	    say $out join "\t", $loc, $source, 'repeat_region', $s, $e, '.', '?', '.', "ID=$rreg";
	    say $out $cop_feats;
	    delete $feature->{$ltr};
	    push @lengths, $len;
	    $pdoms++ if $has_pdoms;
        }
        undef $cop_feats;
        $is_copia = 0;
        $has_pdoms  = 0;
    }
    close $out;

    say STDERR join "\t", "copia_count", "min_length", "max_length", "mean_length", "elements_with_protein_matches";
    my $stat = Statistics::Descriptive::Full->new;
    $stat->add_data(@lengths);
    my $min   = $stat->min;
    my $max   = $stat->max;
    my $mean  = $stat->mean;
    my $count = $stat->count;
    say STDERR join "\t", $count, $min, $max, sprintf("%.2f", $mean), $pdoms;
}

sub write_unclassified_ltrs ($feature, $header, $gff) {
    my @lengths;
    my $unc_feats;
    my $is_unclass = 0;
    my $has_pdoms = 0;
    my $pdoms = 0;
    
    my ($name, $path, $suffix) = fileparse($gff, qr/\.[^.]*/);
    my $outfile = $name."_unclassified.gff3";
    open my $out, '>>', $outfile;
    say $out $header;

    for my $tir (nsort_by { m/repeat_region(\d+)/ and $1 } keys %$feature) {
        my ($rreg, $s, $e) = split /\./, $tir;
        my $len = ($e - $s) + 1;
        my $region = @{$feature->{$tir}}[0];
        my ($loc, $source) = (split /\|\|/, $region)[0..1];
        for my $feat (@{$feature->{$tir}}) {
            my @feats = split /\|\|/, $feat;
            $feats[8] =~ s/\s\;\s/\;/g;
            $feats[8] =~ s/\s+/=/g;
            $feats[8] =~ s/\s+$//;
            $feats[8] =~ s/=$//;
            $feats[8] =~ s/=\;/;/g;
            $feats[8] =~ s/\"//g;
	    $has_pdoms = 1 if $feats[2] =~ /protein_match/;
	    say $out join "\t", $loc, $source, 'repeat_region', $s, $e, '.', '?', '.', "ID=$rreg";
            say $out join "\t", @feats;
        }
	delete $feature->{$tir};
	push @lengths, $len;
	$pdoms++ if $has_pdoms;
	$has_pdoms = 0;
    }
    close $out;

    say STDERR join "\t", "unclassified_count", "min_length", "max_length", "mean_length", "elements_with_protein_matches";
    my $stat = Statistics::Descriptive::Full->new;
    $stat->add_data(@lengths);
    my $min   = $stat->min;
    my $max   = $stat->max;
    my $mean  = $stat->mean;
    my $count = $stat->count;
    say STDERR join "\t", $count, $min, $max, sprintf("%.2f", $mean), $pdoms;
}

sub get_source {
    my ($ref) = @_;

    #dd $ref and exit;
    for my $feat (@$ref) {
	for my $rfeat (@$feat) {
	    my @feats = split /\|\|/, $rfeat;
	    return ($feats[0], $feats[1]);
	}
    }
}
