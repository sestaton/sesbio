#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;
use Bio::DB::Taxonomy;

my $db = Bio::DB::Taxonomy->new(-source => 'entrez');
my $taxonid = $db->get_taxonid('Homo sapiens');
my $taxon = $db->get_taxon(-taxonid => $taxonid);

say "Taxon ID is ", $taxon->id;
say "Scientific name is ", $taxon->scientific_name;
say "Rank is ", $taxon->rank;
say "Division is ", $taxon->division;

if (defined $taxonid) {
    my $node = $db->get_Taxonomy_Node(-taxonid => $taxonid);
    my $kingdom = $node;
    for (1..25) {
	$kingdom = $db->get_Taxonomy_Node(-taxonid => $kingdom->parent_id);
    }
    say "Kingdom is ",$kingdom->scientific_name;
}


