#!/usr/bin/env perl

use 5.014;
use strict;
use warnings;
use LWP::UserAgent;
use HTML::TreeBuilder;
use HTML::TableExtract;
use XML::LibXML;
use Getopt::Long;
use Data::Dump qw(dd dump);

my $usage = "perl $0 -st search_term [-ft] [fs] [-r] [-hmms]\n";
my $search_term;
my $fetch_results;
my $fetch_hmms;
my $family_term;
my $filter_search;

GetOptions(
	   'st|term=s'        => \$search_term,
	   'r|results'        => \$fetch_results,
	   'hmms'             => \$fetch_hmms,
           'ft|family=s'      => \$family_term,
	   'fs|filter_search' => \$filter_search,
	   );

say $usage and exit(1) if !$search_term && !$family_term;
say "\nERROR: Can only choose a search term or a family term, not both." and exit(1) if $search_term and $family_term;

$family_term = ucfirst($family_term) if $family_term;
my $pfam_response = qq{pfam_sanger_search_res.html};
my $ua = LWP::UserAgent->new;
my $tree = HTML::TreeBuilder->new;
$ua->env_proxy;

my ($response, $results);
if ($search_term) {
    $response = $ua->get( "http://pfam.sanger.ac.uk/search/keyword?query=$search_term" );
}
elsif ($family_term) {
    $response = $ua->get( "http://pfam.sanger.ac.uk/family/$family_term?output=xml" );
}

die "error: failed to retrieve XML: " . $response->status_line . "\n"
    unless $response->is_success;

#
# parse the results
#
if ($search_term) {
    $results = parse_search_term($search_term, $response, $pfam_response);
}
elsif ($family_term) {
    $results = parse_family_term($response);
}

#
# print the results
#
if ($fetch_results) {
    dd $results;
}

if ($fetch_hmms) {
    my $to_get = scalar(keys %$results);
    say "========== Will attempt to get $to_get HMMs from the Pfam Database at the Sanger Institute.";
    for my $family (keys %$results) {
	fetch_hmm_files($family);
    }
}

#unlink $pfam_response;
exit;

#
# subroutines
#
sub parse_search_term {
    my ($search_term, $response, $pfam_response) = @_;

    my %results;
    open my $out, '>', $pfam_response or die "\nERROR: Could not open file: $!\n";
    say $out $response->content;
    close $out;
    $tree->parse_file($pfam_response);

    my $te = HTML::TableExtract->new( attribs => { id => q{resultTable} } );
    $te->parse_file($pfam_response);

    for my $ts ($te->tables) {
	for my $row ($ts->rows) {
	    my @elem = grep { defined } @$row;
	    next if $elem[0] =~ /^Original/;
	    if ($filter_search) {
		if ($elem[3] =~ /\Q$search_term\E/) {
		    $results{$elem[1]} = { $elem[2] => $elem[3] };
		}
	    }
	    else {
		$results{$elem[1]} = { $elem[2] => $elem[3] };
	    }
	}
    }
    return \%results;
}

sub parse_family_term {
    my ($response) = @_;

    my %results;
    my $xml = $response->content;

    my $xml_parser = XML::LibXML->new();
    my $dom = $xml_parser->parse_string( $xml );

    my $root = $dom->documentElement();
    my ( $entry ) = $root->getChildrenByTagName( 'entry' );

    my $acc  = $entry->getAttribute('accession');
    my $id   = $entry->getAttribute('id');
    my $desc = $entry->getAttribute('description');
    #print 'accession: ' . $entry->getAttribute( 'accession' ) . "\n";
    say join "\t", $acc, $id, $desc; # for debug
    #$results{$acc} = { $id => $desc };
    #return \%results;
}

sub fetch_hmm_files {
    my ($family) = @_;

    my $ua = LWP::UserAgent->new;
    my $fam_hmm = $family.".hmm";

    # HMM: http://pfam.sanger.ac.uk/family/PF01498/hmm
    my $urlbase = "http://pfam.sanger.ac.uk/family/$family/hmm";
    my $response = $ua->get($urlbase);

    unless ($response->is_success) {
        die "Can't get url $urlbase -- ", $response->status_line;
    }

    say "========== Fetching HMM for $family.";
    open my $out, '>', $fam_hmm or die "\nERROR: Could not open file: $!\n";
    say $out $response->content;
    close $out;
}

