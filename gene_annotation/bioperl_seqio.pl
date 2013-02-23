#!/usr/bin/env perl

use Bio::SeqIO;

my $usage = "$0 infile\n";
my $infile = shift or die $usage;

my $seq_in = Bio::SeqIO->new(-file => $infile, -format => 'fasta');

my $seqct = 0;

while(my $seq = $seq_in->next_seq) {
    $seqct++;
}

print "There are $seqct sequences in $infile.\n";
