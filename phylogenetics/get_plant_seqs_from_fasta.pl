#!/usr/bin/env perl

## Take a fasta file with species names in the header,
## and pull out only the plant species.
##
## NB: This script was used as part of a project to discover centromeric repeats,
##     and it allows you to look at only those derived from plants.
##     The input data to this script was the supplemental fasta file
##     from: http://genomebiology.com/2013/14/1/R10
use 5.010;
use strict;
use warnings;
use Bio::SeqIO;
use Bio::DB::Taxonomy;

my $usage   = "USAGE: $0 seqs.fasta plant_seqs.fasta\n";
my $infile  = shift or die $usage;
my $outfile = shift or die $usage;
my $seqio   = Bio::SeqIO->new(-file => $infile,     -format => 'fasta');
my $seqout  = Bio::SeqIO->new(-file => ">$outfile", -format => 'fasta');
my $db      = Bio::DB::Taxonomy->new(-source => 'entrez');

while (my $seq = $seqio->next_seq) {
    my $id = $seq->id;
    my ($genus, $species) = split /_/, $id, 2;
    my $gen = ucfirst($genus); 
    if ($gen =~/Brachypdium/) { # typo in file (not mine)
	$gen =~s/Brachypdium/Brachypodium/;
    }
    my $taxonid = $db->get_taxonid("$gen $species");
    my $taxon   = $db->get_taxon(-taxonid => $taxonid);
    next unless length($species) > 3; # filter "sp." which is unknown species
    next unless defined $taxon;

    my $div = $taxon->division;
    next unless defined $div;
    if ($div =~ /Plants/i) {
	$seqout->write_seq($seq);
    }
}


