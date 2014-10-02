#!/usr/bin/env perl

## Demonstrate two ways to get taxon information from NCBI
## through BioPerl.

use 5.010;
use strict;
use warnings;

db_taxon();
db_eutils();

sub db_taxon {
    require Bio::DB::Taxonomy;

    my $db      = Bio::DB::Taxonomy->new(-source => 'entrez');
    my $taxonid = $db->get_taxonid('Helianthus annuus');
    my $taxon   = $db->get_taxon(-taxonid => $taxonid);

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
	say "Kingdom is ", $kingdom->scientific_name;
    }
}

sub db_eutils {
    require Bio::DB::EUtilities;
 
    my $id = 4232;
 
    my $factory = Bio::DB::EUtilities->new(-eutil => 'esummary',
					   -email => 'mymail@foo.bar',
					   -db    => 'taxonomy',
					   -id    => $id );
 
    my ($name) = $factory->next_DocSum->get_contents_by_name('ScientificName');
    
    say $name;
}

