#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;
use Bio::SeqIO;
use List::MoreUtils qw(any);

my @ids = (1, 5, 10); 
my $num = 0; 

my $seqio = Bio::SeqIO->new(-fh => \*STDIN); 

while (my $seq = $seqio->next_seq) { 
    $num++; 
    say join "\n", ">".$seq->id, $seq->seq 
	if any { $_ == $num } @ids;
}