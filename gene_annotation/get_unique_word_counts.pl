#!/usr/bin/perl -w

use strict;
#use Data::Dumper;
use Getopt::Long;

my $usage = "USAGE: get_unique_word_counts.pl -i inreport -o outreport

This script takes as input a list of words that may be separated by spaces;
the words are evaluated by unique lines.";

my $infile;
my $outfile;

GetOptions(
           'i|infile=s'     => \$infile,
           'o|outfile=s'    => \$outfile,
           );

# open the infile or die with a usage statement
die $usage if !$infile or !$outfile;

open(my $in, '<', $infile) or die "\nERROR: Could not open file: $infile\n";
open(my $out, '>', $outfile) or die "\nERROR: Could not open file: $outfile\n";

#
# comments must be removed or they will be counted
#
my @words = map +(split "\n")[0], <$in>;

my %seen = ();
my @unique_words = grep { ! $seen{$_} ++ } @words;   # preserves the order of elements
close($in);

my $unique = @unique_words;
my $query = @words;

print "\n","There are: ", $query, " total words.\n";
print "\n","There are: ", $unique, " unique words.\n\n";

count_unique ( @words );

sub count_unique {

    my @array = @_;
    my %count;
    map { $count{$_}++ } @array;

    map {print $out $_."\t".${count{$_}}."\n"} sort keys(%count);

}

close($out);

exit;
