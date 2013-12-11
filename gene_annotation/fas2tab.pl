#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;

my $usage = "$0 infile > out\n";
my $infile = shift or die $usage;

open my $in, '<', $infile or die "\nERROR: Could not open file: $infile\n";

{
    local $/ = '>';

    while (my $line = <$in>) {
        chomp $line;
        my ($seqid, @seqparts) = split /\n/, $line;
        my $seq = join '', @seqparts;
        next unless defined $seqid && defined $seq;
	say join "\t", $seqid, $seq;
    }
}
