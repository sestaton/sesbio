#!/usr/bin/env perl

## This script is used for benchmarking;
## NB: As of this writing, there is memory leak in Bio::Kseq

## Update: memory leak has been fixed: https://github.com/cjfields/Bio-Kseq/issues/1

use strict;
use warnings;
use Bio::Kseq;

my $usage = "$0 infile\n";
my $infile = shift or die $usage;

my $kseq = Bio::Kseq->new($infile);
my $it = $kseq->iterator;
my $ct = 0;

while (my $seq = $it->next_seq) {
    $ct++;
    #print join "\n", ">".$seq->{name}, $seq->{seq},"\n";
}

print "There are $ct sequences in $infile.\n";

