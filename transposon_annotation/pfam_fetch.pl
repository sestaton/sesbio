#!/usr/bin/env perl

=head1 NAME 
                                                                       
pfam_fetch.pl - Fetch data or HMMs for any search term 

=head1 SYNOPSIS    
 
 perl pfam_fetch.pl -st transposase -r -fs

=head1 DESCRIPTION
                                                                   
This is a web client for fetching information, or HMMs, about certain terms 
from the Pfam database at the Sanger Institute. The typical usage would be to
search a general term and print the results. E.g., let's say you wanted to download
the HMM for RVP but did not know the Pfam family name. First, search that term:

    $ perl pfam_fetch.pl -st RVP -r 
    {
      PF00072 => { Response_reg => "Response regulator receiver domain" },
      PF00077 => { RVP => "Retroviral aspartyl protease" },
      PF00478 => { IMPDH => "IMP dehydrogenase / GMP reductase domain" },
      PF03505 => { Clenterotox => "Clostridium enterotoxin" },
      PF08284 => { RVP_2 => "Retroviral aspartyl protease" },
      PF09668 => { Asp_protease => "Aspartyl protease" },
      PF13650 => { Asp_protease_2 => "Aspartyl protease" },
    }

At this point, it is clear that not all of the results are exactly what you want. So, you can
narrow the results by specifying the family you want.

    $ perl pfam_fetch.pl -ft RVP -r
    { PF00077 => "RVP" }

If desired, you can use this information to get the HMM for that model.
    
    $ perl pfam_fetch.pl -ft RVP --hmms
    ========== Will attempt to get 1 HMMs from the Pfam Database at the Sanger Institute.
    ========== Fetching HMM for PF00077.
    ========== Done. Downloaded 1/1 models.

=head1 DEPENDENCIES

This client uses LWP::UserAgent to perform a request, XML::LibXML, 
HTML::TreeBuilder, and HTML::TableExtract to parse the response.

Tested with:

=over 2

=item *

Perl 5.18.0 (on Red Hat Enterprise Linux Server release 5.9 (Tikanga))

=item *

Perl 5.16.0 (on Mac OS X 10.6.8 (Snow Leopard))

=back

=head1 SEE ALSO

HMMER2GO: https://github.com/sestaton/HMMER2GO

This project includes the same functionality as this script and is maintained.

=head1 LICENSE

Copyright (C) 2013-2019 S. Evan Staton

This program is distributed under the MIT (X11) License: http://www.opensource.org/licenses/mit-license.php

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

=head1 AUTHOR 

S. Evan Staton                                                

=head1 CONTACT
 
evan at evanstaton dot com

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

use 5.010;
use strict;
use warnings;
use File::Basename;
use LWP::UserAgent;
use HTML::TreeBuilder;
use HTML::TableExtract;
use XML::LibXML;
use Getopt::Long;
use Pod::Usage;
#use Data::Dump::Color;

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
    say "\nERROR: A search term or Pfam family must be given. Exiting.";
    usage();
    exit(1);
}

if ($search_term && $family_term) {
    say "\nERROR: Choose a search term or a family term, not both. Exiting.";
    usage();
    exit(1);
}


$family_term      = ucfirst($family_term) if $family_term;
my $pfam_response = qq{pfam_sanger_search_res.html};

my $ua   = LWP::UserAgent->new;
my $tree = HTML::TreeBuilder->new;
$ua->env_proxy;

my ($response, $results);
if ($search_term) {
    $response = $ua->get("http://pfam.xfam.org/search/keyword?query=$search_term");
}
elsif ($family_term) {
    $response = $ua->get("http://pfam.xfam.org/family/$family_term?output=xml");
    say $response->content;
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
dd $results if $fetch_results;

my $got_ct = 0;
if ($fetch_hmms) {
    my $to_get = scalar(keys %$results);
    say "========== Will attempt to get $to_get HMMs from the Pfam Database at the Sanger Institute.";
    for my $family (keys %$results) {
	fetch_hmm_files($family);
	$got_ct++
    }
    say "========== Done. Downloaded $got_ct/$to_get models.";
}

unlink $pfam_response;
exit;

#
# methods
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
    my $dom        = $xml_parser->parse_string( $xml );

    my $root      = $dom->documentElement();
    my ( $entry ) = $root->getChildrenByTagName( 'entry' );

    my $acc  = $entry->getAttribute('accession');
    my $id   = $entry->getAttribute('id');

    my @desc_entry = $root->getElementsByTagName( 'description' );
    my $desc = ($desc_entry[0]->childNodes)[1]->nodeValue;
    $desc =~ tr/\n//d;

    ##TODO: Get other information, e.g., go terms
    #my @go_entry = $root->getElementsByTagName( 'go_terms' );
    #my $go1 = $go_entry[0]->nodeType;
    #my $go2 = $go_entry[0]->firstChild->nodeValue;
    #my $go3 = ($go_entry[0]->childNodes)[1]->nodeType;
    #my $go4 = ($go_entry[0]->childNodes)[1]->nodeValue;
    #say join "\n", $go1, $go2, $go3, $go4;

    $results{$acc} = { $id => $desc };

    return \%results;
}

sub fetch_hmm_files {
    my ($family) = @_;

    my $ua       = LWP::UserAgent->new;
    my $fam_hmm  = $family.".hmm";

    my $urlbase  = "http://pfam.sanger.ac.uk/family/$family/hmm";
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

