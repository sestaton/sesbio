#!/usr/bin/env perl

use strict;
use warnings;
use Bio::DB::Taxonomy;

my $db = Bio::DB::Taxonomy->new(-source => 'entrez');
my $taxonid = $db->get_taxonid('Homo sapiens');
my $taxon = $db->get_taxon(-taxonid => $taxonid);

print "Taxon ID is ", $taxon->id, "\n";
print "Scientific name is ", $taxon->scientific_name, "\n";
print "Rank is ", $taxon->rank, "\n";
print "Division is ", $taxon->division, "\n";

if (defined $taxonid) {
    my $node = $db->get_Taxonomy_Node(-taxonid => $taxonid);
    my $kingdom = $node;
    for (1..25) {
	$kingdom = $db->get_Taxonomy_Node(-taxonid => $kingdom->parent_id);
    }
    print "Kingdom is ",$kingdom->scientific_name,"\n";
}


