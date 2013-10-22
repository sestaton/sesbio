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

my $usage = "perl $0 -t search_term [-r] [-hmms]\n";
my $term;
my $fetch_results;
my $fetch_hmms;

GetOptions(
	   't|term=s'   => \$term,
	   'r|results'  => \$fetch_results,
	   'hmms'       => \$fetch_hmms,
	   );

say $usage and exit(1) if !$term;

$term = ucfirst($term);
my $pfam_response = qq{pfam_sanger_search_res.html};
my $ua = LWP::UserAgent->new;
my $tree = HTML::TreeBuilder->new;
$ua->env_proxy;
my %results;

my $res = $ua->get( "http://pfam.sanger.ac.uk/search/keyword?query=$term" );

die "error: failed to retrieve XML: " . $res->status_line . "\n"
    unless $res->is_success;

#
# open and parse the results
#
open my $out, '>', $pfam_response or die "\nERROR: Could not open file: $!\n";
say $out $res->content;
close $out;
$tree->parse_file($pfam_response);

my $te = HTML::TableExtract->new( attribs => { id => q{resultTable} } );
$te->parse_file($pfam_response);

for my $ts ($te->tables) {
    for my $row ($ts->rows) {
        my @elem = grep { defined } @$row;
	next if $elem[0] =~ /^\QOriginal\E/;
	$results{$elem[1]} = { $elem[2] => $elem[3] };
    }
}

if ($fetch_results) {
    dd \%results;
}

if ($fetch_hmms) {
    for my $family (keys %results) {
	fetch_hmm_files($family);
    }
}

unlink $pfam_response;
exit;

#
# subroutines
#
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

    open my $out, '>', $fam_hmm or die "\nERROR: Could not open file: $!\n";
    say $out $response->content;
    close $out;
}

