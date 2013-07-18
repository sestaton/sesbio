#!/usr/local/bin/perl -w
use strict;
use Bio::SeqIO;

if (@ARGV != 2) {    die "USAGE: embl2gb.pl  embl_file gb_file  \n"; }

my $seqio = Bio::SeqIO->new('-format' => 'embl', '-file' => "$ARGV[0]");
my $seqout = new Bio::SeqIO('-format' => 'genbank', '-file' => ">$ARGV[1]");
while( my $seq = $seqio->next_seq) {
  $seqout->write_seq($seq)
}
