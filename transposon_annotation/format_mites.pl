#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;
use Bio::SeqIO;

my $seqio = Bio::SeqIO->new(-fh => \*STDIN);

my $count = 0;

while (my $seq = $seqio->next_seq) {
    $count++;
    my $sequence = $seq->seq;
    $sequence =~ s/.{60}\K/\n/g;
    say join "\n", ">$count"."_".$seq->id."_length=".$seq->length, $sequence;
}

