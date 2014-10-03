#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;
use autodie qw(open);
use File::Basename;
use File::Path qw(make_path);
use File::Spec;
use LWP::UserAgent;
use XML::LibXML;
use HTML::TableExtract;
use Getopt::Long;

my $usage = "$0 -t term term term\n";
my @keywords;
my $outfile;

GetOptions(
	   't|term=s{1,}' => \@keywords,
           'o|outfile=s'  => \$outfile,
	   );

die $usage if !@keywords || !$outfile;

search_by_keyword(\@keywords, $outfile);

sub search_by_keyword {
    my ($keywords, $outfile) = @_;

    my $keyword = join "+", @$keywords;

    my $ua = LWP::UserAgent->new;
    my $urlbase  = "http://pfam.xfam.org/search/keyword?query=$keyword&submit=Submit";
    my $response = $ua->get($urlbase);

    unless ($response->is_success) {
	die "Can't get url $urlbase -- ", $response->status_line;
    }

    my $pfamxml = "pfam_search_$keyword".".xml";
    open my $pfout, '>', $pfamxml;
    say $pfout $response->content;
    close $pfout;

    my ($resultnum, $dbnum) = get_search_results($keyword, $pfamxml);

    if ($resultnum > 1) {
	my $dirname = $keyword; # use expressive variable names

	open my $out, '>', $outfile;
	say "Found $resultnum HMMs for $keyword in $dbnum databases. HMMs can be found in the directory: $dirname.";
	say $out "Accession\tID\tDescription\tPfam";
	make_path($dirname, {verbose => 0, mode => 0771,});
	my $te = HTML::TableExtract->new( headers => [qw(Accession ID Description Pfam Seq_Info Interpro)] );
	$te->parse_file($pfamxml);
	
	for my $ts ($te->tables) {
	    for my $row ($ts->rows) {
		my @elem = grep { defined } @$row;
		#my ($accession, $id, $descripton, $pfam, $seqinfo, $interpro) = @elem;
		fetch_hmm($dirname, \@elem, $out);
	    }
	}
	close $out;
    }
    unlink $pfamxml;
}
    
sub get_search_results {
    my ($keyword, $pfamxml) = @_;
    my ($resultnum, $dbnum);

    $keyword =~ s/\+/ /g;
    open my $in, '<', $pfamxml;
    while (<$in>) {
	if (/We found \<strong\>(\d+)\<\/strong\> unique results/) {
	    $resultnum = $1;
	}
	if (/&quot\;\<em\>$keyword<\/em\>&quot\;\)\, in \<strong\>(\d+)\<\/strong\>/) {
	    $dbnum = $1;
	}
    }
    close $in;

    return $resultnum, $dbnum;
}

sub fetch_hmm {
    my ($dir, $elem, $out) = @_;

    my ($accession, $id, $descripton, $pfam, $seqinfo, $interpro) = @$elem;
    my $ua = LWP::UserAgent->new;
    my $urlbase  = "http://pfam.xfam.org/family/$accession/hmm";
    my $response = $ua->get($urlbase);

    unless ($response->is_success) {
        die "Can't get url $urlbase -- ", $response->status_line;
    }

    say $out join "\t", @$elem[0..3];
    my $hmmfile = File::Spec->catfile($dir, $accession.".hmm");
    open my $hmmout, '>', $hmmfile;
    say $hmmout $response->content;
    close $hmmout;
}
