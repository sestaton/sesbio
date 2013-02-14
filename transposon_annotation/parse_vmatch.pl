#!/usr/bin/env perl

use strict;
use warnings;
#use autodie qw(open);
use feature 'say';

my $ratio; # for taking as an option
my $repeat_ratio = defined($ratio) ? $ratio : "0.80";
my ($match_ct, $simple_match_ct) = (0, 0);

#open(my $in, '<', *STDIN); # or say "ERROR: No input found (This script reads from STDIN, so redirection or a pipe). Exiting." and exit(1);

{
    local $/ = '';
    
    while(<>) {
	$match_ct++;
	chomp;
	next if /^#/;
	my @record = split /\n/;
	my ($match_coords, $subj, $qry) = map { split /\n/ } @record;
	next unless defined $match_coords && defined $subj && defined $qry; # for some reason, vmatch prints only the alignment if it's not a good match 
                                                                            # (that means $qry is not defined but I check all for other corner cases like this)
	$match_coords =~ s/^\s+//;
	my ($slength, $sname, $srpos, $mtype, $qlength, $qname, $qrpos, $dist, $evalue, $score, $pid) = split /\s+/, $match_coords;
	my ($query, $query_string, $qend) = split /\s+/, $qry;
	my ($mono_ratio, $di_ratio) = filter_simple(\$query_string, \$qlength);
	if ($$mono_ratio >= $repeat_ratio || $$di_ratio >= $repeat_ratio) {
	    $simple_match_ct++;
	    #say "$$mono_ratio Monoratio or $$di_ratio Diratio for $query_string is too simple";                
	    say join "\t", $sname, $$mono_ratio, $$di_ratio;
	}         
    }
    my $simple_ratio = sprintf("%.2f",$simple_match_ct/$match_ct);
    say "$simple_ratio of simple k-mers for $match_ct matches";
}
#close($in);

#
# Subs
#
sub filter_simple {
    
    my ($match_string, $merlen) = @_;
    my %di = ('AA' => 0, 'AC' => 0, 'AG' => 0, 'AT' => 0, 'CA' => 0, 'CC' => 0, 'CG' => 0, 'CT' => 0, 'GA' => 0, 'GC' => 0, 'GG' => 0, 'GT' => 0, 'TA' => 0, 'TC' => 0, 'TG' => 0,'TT' => 0);
    my %mono = ('A' => 0, 'C' => 0, 'G' => 0, 'T' => 0);
    my ($mono_ct, $di_ct, $mono_ratio, $di_ratio) = (0, 0, 0, 0);

    #my %simpleseqs;

    for my $mononuc (keys %mono) {
	while ($$match_string =~ /$mononuc/ig) { $mono_ct++; }
        $mono_ratio = sprintf("%.2f",$mono_ct/$$merlen);
 	#if ($monoratio >= $repeat_ratio) {
	#    $simpleseqs{$seq} = $monoratio;
	#}
	$mono_ct = 0;
    }
    
    for my $dinuc (keys %di) {
	while ($$match_string =~ /$dinuc/ig) { $di_ct++; }
	$di_ratio = sprintf("%.2f",$di_ct*2/$$merlen);
	#if ($diratio >= $repeat_ratio) {
	#    $simpleseqs{$seq} = $diratio;
	#}
	$di_ct = 0;
    }

    return(\$mono_ratio, \$di_ratio);

}









