#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;
use LWP::UserAgent;
use HTML::TableExtract;
use Data::Dump::Color;

my $outfile  = 'tair_tetable.html';
my $ua       = LWP::UserAgent->new;
my $urlbase  = 'https://www.arabidopsis.org/servlets/processor?type=transposonfamily&update_action=browse';
my $response = $ua->get($urlbase);

unless ($response->is_success) {
    die "Can't get url $urlbase -- ", $response->status_line;
}

open my $out, '>', $outfile or die "\nERROR: Could not open file: $!\n";
say $out $response->content;
close $out;

say join "\t", "Family", "Superfamily", "Count";

my $te = HTML::TableExtract->new( attribs => { border => 0 } );
$te->parse_file($outfile);

for my $ts ($te->tables) {
    for my $row ($ts->rows) {
	my @elem = grep { defined } @$row;
	for my $str (@elem) {
	    $str =~ s/^\s+//u;
	}
	@elem = grep { /\S+/ } @elem;
	next if $elem[0] =~ /^.{1,1}$|^Family|back to top/i;
	say join "\t", @elem;
    }
}

unlink $outfile;
