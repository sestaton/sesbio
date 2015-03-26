#!/usr/bin/env perl

## A pure-Perl parser for eSearch results. This is easier/faster
## to install than libxml2 (required by XML::LibXML) and it is fast enough.

use 5.010;
use strict;
use warnings;
use XML::Twig;

my $xml = do { local $/; <DATA> };

my $parser = XML::Twig->new;
$parser->parse($xml);

my @nodes = $parser->findnodes( '/eSearchResult/IdList/Id' );
my $id = pop(@nodes)->text();
say $id;

__DATA__
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE eSearchResult PUBLIC "-//NLM//DTD esearch 20060628//EN" "http://eutils.ncbi.nlm.nih.gov/eutils/dtd/20060628/esearch.dtd">
<eSearchResult><Count>1</Count><RetMax>1</RetMax><RetStart>0</RetStart><IdList>
<Id>4232</Id>
</IdList><TranslationSet/><TranslationStack>   <TermSet>    <Term>Helianthus annuus[All Names]</Term>    <Field>All Names</Field>    <Count>1</Count>    <Explode>N</Explode>   </TermSet>   <OP>GROUP</OP>  </TranslationStack><QueryTranslation>Helianthus annuus[All Names]</QueryTranslation></eSearchResult>
