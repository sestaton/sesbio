#!/usr/bin/env perl

use Bio::Kseq;

my $usage = "$0 infile\n";
my $infile = shift or die $usage;

my $kseq = Bio::Kseq->new($infile);
my $it = $kseq->iterator;
my $ct = 0;
#my $sample_seq;
while (my $seq = $it->next_seq) {
    $ct++;
    #$sample_seq = $seq; # always grab the last seq
    #print join("\n",(">".$seq->{name},$seq->{seq})),"\n";
}

print "There are $ct sequences in $infile.\n";

