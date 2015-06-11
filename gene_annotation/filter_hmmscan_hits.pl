#!/usr/bin/env perl

## TODO: use hash slice for splitting report and list-moreutils 'any' for
## finding matches

use 5.010;
use strict;
use warnings;
use Getopt::Long;

my $usage = "\n$0 -i list -j hmmscan_report -o outfile\n";
my $infile;
my $outfile;
my $hmmfile;

GetOptions(
           'i|infile=s'    => \$infile,
           'j|hmmfile=s'   => \$hmmfile,
           'o|outfile=s'   => \$outfile,
           );

if (!$infile) {
    die "\nERROR: no infile found.\n",$usage;
}
if (!$hmmfile) {
    die "\nERROR: No hmmscan file found.\n",$usage;
}
if (!$outfile) {
    die "\nERROR: No outfile found.\n",$usage;
}

open my $in, '<', $infile or die "\nERROR: Could not open file: $infile\n";
open my $hmm, '<', $hmmfile or die "\nERROR: Could not open file: $hmmfile\n";

open my $out, '>', $outfile or die "\nERROR: Could not open file: $outfile\n";

my %ids;
while (<$in>) {
    chomp;
    next if /^\#/;
    $ids{$_} = 1;
}
close $in;

while (<$hmm>) {
    chomp;
    next if /^\#/;
    my ($target_name, $accession, $query_name, $accession_q, $E_value_full, $score_full, 
	$bias_full, $E_value_best, $score_best, $bias_best, $exp, $reg, $clu, $ov, $env, 
	$dom,$rev, $inc, $description_of_target) = split /\t/;

    for my $key (sort keys %ids) {
	if ($key eq $query_name) {
	    say $out join "\t", $target_name, $accession, $query_name, $accession_q, $E_value_full, 
	        $score_full, $bias_full, $E_value_best, $score_best, $bias_best, $exp, 
	        $reg, $clu, $ov, $env, $dom, $rev, $inc, $description_of_target;
	}
    }
}
close $hmm;
close $out;
