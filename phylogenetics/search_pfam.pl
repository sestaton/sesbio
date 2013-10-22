#!/usr/bin/env perl

use strict;
use warnings;
use LWP::UserAgent;
use XML::LibXML;

my $usage = "perl $0 family_name\n";
my $family = shift or die $usage;
$family = ucfirst($family);

#print $family and exit;

my $ua = LWP::UserAgent->new;
$ua->env_proxy;

my $res = $ua->get( 'http://pfam.sanger.ac.uk/family/$family?output=xml' );

die "error: failed to retrieve XML: " . $res->status_line . "\n"
    unless $res->is_success;

my $xml = $res->content;

my $xml_parser = XML::LibXML->new();
my $dom = $xml_parser->parse_string( $xml );

my $root = $dom->documentElement();
my ( $entry ) = $root->getChildrenByTagName( 'entry' );

print 'accession: ' . $entry->getAttribute( 'accession' ) . "\n";
