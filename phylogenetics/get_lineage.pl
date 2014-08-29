#!/usr/bin/env perl

##NB:   This demonstrates how to get a taxonomic lineage w/o the use of BioPerl.
##TODO: Build lineage info based on....(the number of entries?)
  
use 5.010;
use strict;
use warnings;
use LWP::Simple;
use XML::LibXML;

my $id    = 4232; # Helianthus annuus
my $esumm = "esumm_$id.xml"; 
 
my $ua = LWP::UserAgent->new;
my $urlbase  = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id=$id";
my $response = $ua->get($urlbase);

unless ($response->is_success) {
    die "Can't get url $urlbase -- ", $response->status_line;
}

open my $out, '>', $esumm or die "\nERROR: Could not open file: $!\n";
say $out $response->content;
close $out;

my $parser = XML::LibXML->new;
my $doc    = $parser->parse_file($esumm);

for my $node ( $doc->findnodes('//TaxaSet/Taxon') ) {
    my ($lineage) = $node->findvalue('Lineage/text()');
    my ($family)  = map  { s/\;$//; $_; } 
	            grep { /(\w+aceae)/ } 
                    map  { split /\s+/  } $lineage;
    say $family;
    say $lineage;
}

unlink $esumm;
