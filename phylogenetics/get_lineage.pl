#!/usr/bin/env perl

##TODO: Build lineage info based on....(the number of entries?)
  
use 5.010;
use strict;
use warnings;
use LWP::Simple;
use XML::LibXML;

my $id      = 4232; # Helianthus annuus
my $genus   = 'Helianthus';
my $species = 'annuus';

search_by_name($genus, $species);

#
# methods
#
sub search_by_name {
    my ($genus, $species) = @_;

    my $id = _fetch_id_for_name($genus, $species);
    say join "\t", $genus, $species, $id;
    
    _get_lineage_for_id($id);
}

sub _get_lineage_for_id {
    my ($id) = @_;
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
	say "Family: $family";
	say "Full taxonomic lineage: $lineage";
    }
    
    unlink $esumm;
}

sub _fetch_id_for_name {
    my ($genus, $species) = @_;

    my $esearch = "esearch_$genus"."_"."$species.xml";
    my $ua = LWP::UserAgent->new;
    my $urlbase  = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=taxonomy&term=$genus%20$species";
    my $response = $ua->get($urlbase);

    unless ($response->is_success) {
	die "Can't get url $urlbase -- ", $response->status_line;
    }

    open my $out, '>', $esearch or die "\nERROR: Could not open file: $!\n";
    say $out $response->content;
    close $out;

    my $id;
    my $parser = XML::LibXML->new;
    my $doc    = $parser->parse_file($esearch);
    
    for my $node ( $doc->findnodes('//eSearchResult/IdList') ) {
	($id) = $node->findvalue('Id/text()');
    }
    
    unlink $esearch;
    
    return $id;
}
