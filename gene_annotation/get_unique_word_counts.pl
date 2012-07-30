#!/usr/bin/perl -w

use strict;
#use Data::Dumper;
use Getopt::Long;

my $usage = "\nUSAGE: $0 -i in -o out <--stats>

This script takes as input a list of words that may be separated by spaces;
the words are evaluated by unique lines

Will also read from a Unix pipe.\n\n";

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

if ($help) {
    print $usage;
    exit(0);
}

# open the infile 
open(my $words, '<', $infile) or die "\nERROR: Cannot open file: $infile\n" if $infile;
open(my $unique_words, '>', $outfile) or die "\nERROR: Cannot open file: $outfile\n" if $outfile;

#
# comments must be removed or they will be counted
#

@names = map +(split "\n")[0], <$words> if $infile;
@names = map +(split "\n")[0], <STDIN> if !$infile;

my %seen = ();
my @unique_names = grep { ! $seen{$_} ++ } @names;   # preserves the order of elements
close($words) if $infile;

my $unique = @unique_names;
my $total = @names;

if ($statistics) {
    print "\nThere are: ", $total, " total words.\n";
    print "\nThere are: ", $unique, " unique words.\n\n";
}

count_unique ( @names );

exit;

#
# subs
#
sub count_unique {

    my @array = @_;
    my %count;
    map { $count{$_}++ } @array;

    if ($outfile) {
	map {print $unique_words $_."\t".${count{$_}}."\n"} sort keys(%count);
	close($unique_words);
    } else {
	map {print $_."\t".${count{$_}}."\n"} sort keys(%count);
    }

}

