#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;
use Getopt::Long;

my $usage = "\nUSAGE: $0 -i in -o out <--stats>

This script takes as input a list of words that may be separated by spaces;
the words are evaluated by unique lines.

Input may be a file or from stdin.\n";

my $infile;
my $outfile;
my $statistics;
my @names;
my $help;

GetOptions(
           'i|infile=s'    => \$infile,
           'o|outfile=s'   => \$outfile,
           'stats'         => \$statistics,
           'h|help'        => \$help,
           );

say $usage and exit(0) if $help;

open my $words, '<', $infile or die "\nERROR: Cannot open file: $infile\n" if $infile;
open my $unique_words, '>', $outfile or die "\nERROR: Cannot open file: $outfile\n" if $outfile;

#
# comments must be removed or they will be counted
#
@names = map +(split "\n")[0], <$words> if $infile;
@names = map +(split "\n")[0], <STDIN> if !$infile;

my %seen = ();
my @unique_names = grep { ! $seen{$_} ++ } @names;   # preserves the order of elements
close $words if $infile;

my $unique = @unique_names;
my $total  = @names;

if ($statistics) {
    say "\nThere are: ", $total, " total words.";
    say "\nThere are: ", $unique, " unique words.\n";
}

count_unique(@names);

exit;
#
# methods
#
sub count_unique {
    my @array = @_;
    my %count;
    map { $count{$_}++ } @array;

    if ($outfile) {
	map { say $unique_words join "\t", $_, ${count{$_}} } sort keys %count;
	close $unique_words;
    } 
    else {
	map { say join "\t", $_, ${count{$_}} } sort keys %count;
    }
}

