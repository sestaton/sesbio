#!/usr/bin/env perl

use v5.10;
use strict; 
use warnings;
use Getopt::Long;
#use Data::Dumper;

my $usage = "\n$0 -i annot -p pfam2go -o outfile <--map>\n";
my $infile;
my $pfam2go; 
my $outfile;
my $mapping;
my $mapfile;
my $map_fh;

GetOptions(
	   'i|infile=s'  => \$infile,
	   'p|pfam2go=s' => \$pfam2go,
	   'o|outfile=s' => \$outfile,
	   'map'         => \$mapping,
	   );

if (!$infile) {
    die "\nERROR: no infile found.\n",$usage;
}
if (!$pfam2go) {
    die "\nERROR: No pfam2go db found.\n",$usage;
}
if (!$outfile) {
    die "\nERROR: No outfile found.\n",$usage;
}
open my $in, '<', $infile or die "\nERROR: Could not open file: $infile\n";
open my $pfams, '<', $pfam2go or die "\nERROR: Could not open file: $pfam2go\n";

open my $out, '>', $outfile or die "\nERROR: Could not open file: $outfile\n";

if ($mapping) {
    $mapfile = $outfile;
    $mapfile =~ s/\..*//g;
    $mapfile .= "_GOterm_mapping.txt";
    open $map_fh, '>', $mapfile or die "\nERROR: Could not open file: $mapfile\n";
}

my %pfamids;
while(<$in>) {
    chomp;
    next if /^\#/;
    my ($target_name, $accession, $query_name, $accession_q, $E_value_full, 
	$score_full, $bias_full, $E_value_best, $score_best, $bias_best, 
	$exp, $reg, $clu, $ov, $env, $dom, $rev, $inc, $description_of_target) = split;
    my $query_eval = join ",", $query_name, $E_value_full, $description_of_target;
    $accession =~ s/\..*//;
    $pfamids{$query_eval} = $accession;
}
close $in;

#print Dumper %pfamids;

#my %mappedterms;
my %goterms;
my $go_ct = 0;
my $map_ct = 0;

while(my $mapping = <$pfams>) {
    chomp $mapping;
    next if $mapping =~ /^!/;
    if ($mapping =~ /Pfam:(\S+) (\S+ \> )(GO\:\S+.*\;) (GO\:\d+)/) {
	my $pf = $1;
	my $pf_name = $2;
	my $pf_desc = $3;
	my $go_term = $4;
	$pf_name =~ s/\s.*//;
	$pf_desc =~ s/\s\;//;
	for my $key (keys %pfamids) { 
	    my ($query, $eval, $desc) = split /\,/, $key;
	    if ($pfamids{$key} eq $pf) {
		say $out join "\t", $query, $pf, $pf_name, $pf_desc, $go_term, $desc;
		if ($mapping) {
		    if (exists $goterms{$query}) {
			$go_ct++ if defined($go_term);
			$goterms{$query} .= ",".$go_term;
		    } else {
			$goterms{$query} = $go_term;
		    }
		}
		last;
	    }
	}
    }
}
close $pfams;
close $out);

if ($mapping) {
    while(my ($key, $value) = each(%goterms)) {
	$map_ct++;
	say $map_fh join "\t", $key, $value;
    }
    say "\n$map_ct query sequences with $go_ct GO terms mapped in file $mapfile.\n";
}






