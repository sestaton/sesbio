#!/usr/bin/perl -w

use strict;
use Bio::SeqIO;

my $usage = "$0 infile > out\n";
my $infile = shift or die $usage;

my $seq_in = Bio::SeqIO->new(-file => $infile,
			     -format => 'fasta');

while(my $seq = $seq_in->next_seq) {
    print join("\t",($seq->id,$seq->seq)),"\n";
}

exit;
