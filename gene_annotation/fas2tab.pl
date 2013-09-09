#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;
use Bio::SeqIO;

my $usage = "$0 infile > out\n";
my $infile = shift or die $usage;

my $seq_in = Bio::SeqIO->new(-file => $infile,
			     -format => 'fasta');

while(my $seq = $seq_in->next_seq) {
    say join "\t", $seq->id, $seq->seq;
}
