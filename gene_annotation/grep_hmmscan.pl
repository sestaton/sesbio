#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;

my $usage = "\n$0 -i list -j hmmscan_report -o outfile\n";
my $infile;
my $outfile;
my $blast;

GetOptions(
           'i|infile=s'  => \$infile,
           'j|blastfile=s' => \$blast,
           'o|outfile=s' => \$outfile,
           );

if (!$infile) {
    die "\nERROR: no infile found.\n",$usage;
}
if (!$blast) {
    die "\nERROR: No hmmscan file found.\n",$usage;
}
if (!$outfile) {
    die "\nERROR: No outfile found.\n",$usage;
}
open my $in, '<', $infile or die "\nERROR: Could not open file: $infile\n";
open my $bl, '<', $blast or die "\nERROR: Could not open file: $blast\n";

open my $out, '>', $outfile or die "\nERROR: Could not open file: $outfile\n";

my %ids;
while(<$in>) {
    chomp;
    next if /^\#/;
    $ids{$_} = 1;
}
close $in;

#my %pfamids;

while(<$bl>) {
    chomp;
    next if /^\#/;
    my ($target_name, $accession, $query_name, $accession_q, $E_value_full, $score_full, 
	$bias_full, $E_value_best, $score_best, $bias_best, $exp, $reg, $clu, $ov, $env, 
	$dom,$rev, $inc, $description_of_target) = split;
    #my $query_eval = join(",",$query_name,$E_value_full,$description_of_target);
    #$accession =~ s/\..*//;
    #$pfamids{$query_eval} = $accession;
    for my $key (sort keys %ids) {
	if ($key eq $query_name) {
	    say $out join "\t", $target_name, $accession, $query_name, $accession_q, 
	    $E_value_full, $score_full, $bias_full, $E_value_best, $score_best, 
	    $bias_best, $exp, $reg, $clu, $ov, $env, $dom, $rev, $inc, $description_of_target;
	}
    }
}
close $bl;
close $out;
