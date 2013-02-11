#!/usr/bin/env perl

use strict;
use warnings;
use Bio::SeqIO;

my $usage = "$0 in out\n";
my $in = shift or die $usage;
my $out = shift or die $usage;

my $seq_in = Bio::SeqIO->new(-file => $in, -format => 'fasta');
my $seq_out = Bio::SeqIO->new(-file => ">$out", -format => 'fasta');

my %seqs;
while (my $seq = $seq_in->next_seq) {
    unless (exists $seqs{$seq->id}) {
	$seq_out->write_seq($seq);
    }
    $seqs{$seq->id} = 1;
}
    
