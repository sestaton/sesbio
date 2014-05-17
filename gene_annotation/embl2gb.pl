#!/usr/bin/env perl

use strict;
use warnings;
use Bio::SeqIO;

my $usage     = "USAGE: embl2gb.pl embl_file gb_file \n";
my $embl_file = shift or die $usage;
my $gb_file   = shift or die $usage;

my $seqio  = Bio::SeqIO->new(-format => 'embl',    -file => $embl_file);
my $seqout = Bio::SeqIO->new(-format => 'genbank', -file => ">$gb_file");

while (my $seq = $seqio->next_seq) {
  $seqout->write_seq($seq)
}
