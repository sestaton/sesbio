#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;
use Bio::Tools::GFF;
use Data::Dump::Color;
use File::Find;
use Cwd;
use Sort::Naturally;
use List::UtilsBy qw(nsort_by);

my $gff = shift;
my $cwd = getcwd();
my @files;
find( sub { push @files, $File::Find::name if -f and /tidy_Ha\d+\.gff3$/ }, $cwd );

my ($expct, $pepct, $mrnact, $cdsct, $genect, $headct) = (0, 0, 0, 0, 0, 1);

for my $file (nsort @files) {
    my ($header, $features) = collect_gff_features($gff, $file);
    say $header if $headct;
    for my $id (nsort_by { m/\w+\.(\d+)\.\d+/ and $1 } keys %$features) {
	my ($parentid, $start, $stop) = split /\./, $id;
	$genect++ if $parentid =~ /gene/;
	$mrnact++ if $parentid =~ /mRNA/;
	$cdsct++  if $parentid =~ /CDS/;
	$pepct++  if $parentid =~ /protein_match/;
	$expct++  if $parentid =~ /expressed_sequence_match/;
	for my $parent (keys %{$features->{$id}}) {
	    my @parent_feats = split /\|\|/, $parent;
	    $parent_feats[8] = _format_attribute($parent_feats[8], $genect, $mrnact, $cdsct, $pepct, $expct);
	    say join "\t", @parent_feats;
	    for my $feat (@{$features->{$id}{$parent}}) {
		my @part_feats =  split /\|\|/, $feat;
		$part_feats[8] = _format_attribute($part_feats[8], $genect, $mrnact, $cdsct, $pepct, $expct);
		say join "\t", @part_feats;
	    }
	}
    }
    $headct = 0;
}

sub collect_gff_features {
    my ($gff, $file) = @_;

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

    my $gffio = Bio::Tools::GFF->new( -file => $file, -gff_version => 3 );

    my ($start, $end, $region, $parent, %features);
  FEATURE:
    while (my $feature = $gffio->next_feature()) {
	if ($feature->primary_tag =~ /protein_match|expressed_sequence_match|gene/) {
	    my @string = split /\t/, $feature->gff_string;
	    ($region) = ($string[8] =~ /ID=?\s+?(protein_match\d+|expressed_sequence_match\d+|gene\d+)/);
	    ($start, $end) = ($feature->start, $feature->end);
	    $parent = join "||", @string;
	}
	next FEATURE unless defined $start && defined $end;
	if ($feature->primary_tag !~ /protein_match|expressed_sequence_match|gene/) {
	    if ($feature->start >= $start && $feature->end <= $end) {
		push @{$features{$region.".".$start.".".$end}{$parent}}, 
		       join "||", split /\t/, $feature->gff_string;
	    }
	}
    }

    return ($header, \%features);
}

sub _format_attribute {
    my ($str, $genect, $mrnact, $cdsct, $pepct, $expct) = @_;

    $str =~ s/\s\;\s/\;/g;
    $str =~ s/\s+/=/g;
    $str =~ s/\s+$//;
    $str =~ s/=$//;
    $str =~ s/=\;/;/g;
    $str =~ s/\"//g;
    
    $str =~ s/gene\d+/gene$genect/ 
	if $str =~ /gene|mRNA/;
    $str =~ s/mRNA\d+/mRNA$mrnact/ 
	if $str =~ /mRNA|exon|three_prime_UTR|five_prime_UTR|CDS/;
    $str =~ s/CDS\d+/CDS$cdsct/ 
	if $str =~ /CDS/;
    $str =~ s/protein_match\d+/protein_match$pepct/ 
	if $str =~ /protein_match/;
    $str =~ s/expressed_sequence_match\d+/expressed_sequence_match$expct/ 
	if $str =~ /expressed_sequence_match/;

    return $str;
}
