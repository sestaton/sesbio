#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;

my $usage = "$0 infile > out OR $0 infile | some_other_program ...\n";
my $infile = shift or die $usage;

open my $in, '<', $infile or die "\nERROR: Could not open file: $infile\n";
    
while (my $line = <$in>) {
    chomp $line;
    my ($header, $seq) = split /\t/, $line;
    #$seq =~ s/(.{60})/$1\n/gs;  
    $seq =~ s/.{60}\K/\n/g;
    say join "\n", ">".$header, $seq;
}

close $in;
