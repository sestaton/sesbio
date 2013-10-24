#!/usr/bin/env perl

=head1 NAME 
                                                                       
pfam_fetch.pl - Fetch data or HMMs for any search term 

=head1 SYNOPSIS    
 
pfam_fetch.pl -st transposase -r -fs

=head1 DESCRIPTION
                                                                   

(...)

=head1 DEPENDENCIES

This client uses LWP::UserAgent to perform a request, XML::LibXML, 
HTML::TreeBuilder, and HTML::TableExtract to parse the response.

Tested with:

(to be added later)

=head1 LICENSE

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

=head1 AUTHOR 

S. Evan Staton                                                

=head1 CONTACT
 
statonse at gmail dot com

=head1 REQUIRED ARGUMENTS

=over 2

=item -st, --search_term

A term to search Pfam for related entries.

=item -ft, --family_term

A Pfam family term to search for a specific entry.

=back

=head1 OPTIONS

=over 2

=item -r, --results

Print a table of search results for introspection.

=item --hmms

Download HMMs for each search result. It is advisable to print the results
first and make sure there are no spurrious results. Also, it may be helpful
to filter the search results (see below).

=item -f, --filter_search

Use the search term to filter the descriptions in the results. This is not always
necessary, but some terms may return dozens of results that are are unannotated.

=item -h, --help

Print a usage statement. 

=item -m, --man

Print the full documentation.

=cut      

use 5.014;
use strict;
use warnings;
use File::Basename;
use LWP::UserAgent;
use HTML::TreeBuilder;
use HTML::TableExtract;
use XML::LibXML;
use Getopt::Long;
use Pod::Usage;
use Data::Dump qw(dd dump);

my $search_term;
my $fetch_results;
my $fetch_hmms;
my $family_term;
my $filter_search;
my $help;
my $man;

GetOptions(
	   'st|term=s'        => \$search_term,
	   'r|results'        => \$fetch_results,
	   'hmms'             => \$fetch_hmms,
           'ft|family=s'      => \$family_term,
	   'fs|filter_search' => \$filter_search,
	   'h|help'           => \$help,
	   'm|man'            => \$man,
	   );

pod2usage( -verbose => 1 ) if $help;
pod2usage( -verbose => 2 ) if $man;

if (!$search_term && !$family_term) {
    say "\nERROR: A search term or Pfam family must be given. Exiting.\n";
    usage();
    exit(1);
}

if ($search_term && $family_term) {
    say "\nERROR: Can only choose a search term or a family term, not both. Exiting.\n";
    usage();
    exit(1);
}

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

unlink $pfam_response;
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
    $results{$acc} = $id;
    #$results{$acc} = { $id => $desc }; # can't get description from family search form
    return \%results;
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

sub usage {
    my $script = basename($0);
  print STDERR <<END

USAGE: $script -st term -r [-ft] [-fs] [--hmms] [-h] [-m]

Required:
    -st|term           :    The name of a domain or protein to search Pfam for similar entries.    

Options:
    -ft|family         :    A specific Pfam family to search for in the Pfam database.
    -fs|filter_search  :    Filter the descriptions of individual results to remove non-specific results.
    -r|results         :    Print a table of results for inspecting different search terms.
    --hmms             :    Download the HMMs for each Pfam family returned by the search.
    -h|help            :    Print a usage message and exit.
    -m|man             :    Print the full documation.

END
}

