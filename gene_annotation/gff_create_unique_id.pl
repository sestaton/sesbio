#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;
use Bio::Tools::GFF;
use File::Find;
use Cwd;
use Sort::Naturally;
use List::UtilsBy qw(nsort_by);

my $gff = shift;
my $cwd = getcwd();
my @files;
find( sub { push @files, $File::Find::name if -f and /tidy_Ha\d+\.gff3$/ }, $cwd );

my ($expct, $pepct, $mrnact, $cdsct, $genect, $threeprutrct, $fiveprutrct, 
    $headct, $ingene) = (0, 0, 0, 0, 0, 0, 0, 1, 1);
my ($has_mrna, $has_cds) = (0, 0);

for my $file (nsort @files) {
    my ($header, $features) = collect_gff_features($gff, $file);
    say $header if $headct;
    for my $id (nsort_by { m/\w+\.(\d+)\.\d+/ and $1 } keys %$features) {
	my ($parentid, $start, $stop) = split /\./, $id;
	$genect++ if $parentid =~ /gene/;
	$pepct++  if $parentid =~ /protein_match/;
	$expct++  if $parentid =~ /expressed_sequence_match/;
	for my $parent (keys %{$features->{$id}}) {
	    my @parent_feats = split /\|\|/, $parent;
	    $parent_feats[8] = 
		_format_parent_attribute($parent_feats[8], $genect, $pepct, $expct);
	    say join "\t", @parent_feats;
	    
	    ($mrnact, $cdsct, $threeprutrct, $fiveprutrct) 
		= _check_part_features(\@{$features->{$id}{$parent}}, $mrnact, $cdsct, $threeprutrct, $fiveprutrct);
	    for my $feat (@{$features->{$id}{$parent}}) {
		my @part_feats = split /\|\|/, $feat;
		$part_feats[8] = 
		    _format_part_attribute($part_feats[2], $part_feats[8], $genect, $mrnact, $cdsct, $threeprutrct, $fiveprutrct);
		say join "\t", @part_feats;
	    }
	}
	say "###";
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

sub _check_part_features {
    my ($feats, $mrnact, $cdsct, $threeprutrct, $fiveprutrct) = @_;

    my ($has_mrna, $has_cds, $has_thrprutr, $has_fiveprutr) = (0, 0);
    for my $feat (@$feats) {
	my @part_feats = split /\|\|/, $feat;
	#$has_mrna = 1 if $part_feats[2] =~ /mRNA/;
	$mrnact++ if $part_feats[2] =~ /mRNA/;
	$has_cds = 1 
	    if $part_feats[2] =~ /CDS/ && $part_feats[8] =~ /ID=?\s+?CDS/;
	$has_thrprutr = 1 
	    if $part_feats[2] =~ /three_prime_UTR/ && $part_feats[8] =~ /ID=?\s+?three_prime_UTR/;
	$has_fiveprutr = 1 
	    if $part_feats[2] =~ /five_prime_UTR/ && $part_feats[8] =~ /ID=?\s+?five_prime_UTR/;
    }
    #$mrnact++ if $has_mrna;
    $cdsct++ if $has_cds;
    $threeprutrct++ if $has_thrprutr;
    $fiveprutrct++  if $has_fiveprutr;

    return ($mrnact, $cdsct, $threeprutrct, $fiveprutrct);
}

sub _format_parent_attribute {
    my ($str, $genect, $pepct, $expct) = @_;

    $str =~ s/\s\;\s/\;/g;
    $str =~ s/\s+/=/g;
    $str =~ s/\s+$//;
    $str =~ s/=$//;
    $str =~ s/=\;/;/g;
    $str =~ s/\"//g;
    
    $str =~ s/gene\d+/gene$genect/ 
	if $str =~ /gene/;
    $str =~ s/protein_match\d+/protein_match$pepct/ 
	if $str =~ /protein_match/;
    $str =~ s/expressed_sequence_match\d+/expressed_sequence_match$expct/ 
	if $str =~ /expressed_sequence_match/;

    return $str;
}

sub _format_part_attribute {
    my ($tag, $str, $genect, $mrnact, $cdsct, $threeprutrct, $fiveprutrct) = @_;

    $str =~ s/\s\;\s/\;/g;
    #$str =~ s/\s+/=/g;
    $str =~ s/\s+$//;
    $str =~ s/=$//;
    $str =~ s/=\;/;/g;
    $str =~ s/\"//g;
    
    $str =~ s/gene\d+/gene$genect/ 
        if $tag =~ /mRNA/;
    $str =~ s/mRNA\d+/mRNA$mrnact/ 
        if $tag =~ /mRNA|exon|three_prime_UTR|five_prime_UTR|CDS/;
    $str =~ s/CDS\d+/CDS$cdsct/ 
        if $tag =~ /CDS/;
    $str =~ s/three_prime_UTR\d+/three_prime_UTR$threeprutrct/
        if $tag =~ /three_prime_UTR/;
    $str =~ s/three_prime_UTR\d+/five_prime_UTR$fiveprutrct/
	if $tag =~ /five_prime_UTR/;

    return $str;
}
