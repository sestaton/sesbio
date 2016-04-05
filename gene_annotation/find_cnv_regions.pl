#!/usr/bin/env perl

use 5.022;
use warnings;
use Set::IntervalTree;
use Data::Dump::Color;
use experimental 'signatures';

my $usage = "$0 cnvs";
my $cnvs  = shift or die $usage;
my $map   = shift or die $usage;

my $lg4_tree  = Set::IntervalTree->new;
my $lg16_tree = Set::IntervalTree->new;

#Ha16 150857390 150861558 LG16_scf7180038271797
#Ha4  174602340 174608711 LG4_scaffold_500
$lg4_tree->insert(  'Ha4',  174602340, 174608711);
$lg16_tree->insert( 'Ha16', 150857390, 150861558);

my %gene_coords = (
    'Ha4'  => $lg4_tree,
    'Ha16' => $lg16_tree,    
);

my $samples = get_samplenames($map);
#dd $samples and exit;
parse_calls($cnvs, $samples, \%gene_coords);

sub get_samplenames ($samplemap) {
    my %samples;

    open my $in, '<', $samplemap;
    my $head = <$in>;

    while (my $line = <$in>) {
        chomp $line;
        my ($resistance, $category, $indiv, $deathdate, $intermating) = split /\t/, $line;
        #Least_ResistantBottom_10%IA1A-45-AugNANA
        $samples{$indiv} = join "||", $resistance, $category, $deathdate, $intermating;
    }
    close $in;

    return \%samples;
}

sub parse_calls ($calls, $samples, $gene_coords) {
    open my $in, '<', $calls;
    while (my $line = <$in>) {
        chomp $line;
	next if $line =~ /^seqnames/;
        my ($index, $chr, $start, $end, $width, $strand, $sample, $median, $mean, $cn) = split /\t/, $line;
        #my ($ref, $s, $e) = split /[:-]/, $f[1];
        my $res = $gene_coords->{$chr}->fetch($start, $end);
        if (@$res) {
	    my $id = $sample;
	    $id =~ s/.*HT\d+_//;
	    $id =~ s/_sort.*//;
	    if (exists $samples->{$id}) {
		#say STDERR "debug: $id";
		my ($type, $what, $na, $date) = split /\|\|/, $samples->{$id};
		#say "CNV on target region: ",@$res;
		say join "\t", $index, $type, $what, $chr, $start, $end, $width, $strand, $sample, $median, $mean, $cn;
	    }
	}
    }
    close $in;
}
