#!/usr/bin/env perl

use strict;
use warnings;
use Bio::SeqIO;
use Getopt::Long;

my $usage = "perl $0 -i seqs_anno.embl -o seqs_anno.gb\n";
my $infile;
my $outfile;

GetOptions(
           'i|infile=s'  => \$infile,
           'o|outfile=s' => \$outfile,
           );

die $usage if !$infile or !$outfile;

my $seqin  = Bio::SeqIO->new(-file => $infile,     -format => 'embl');
my $seqout = Bio::SeqIO->new(-file => ">$outfile", -format => 'genbank');

while (my $seq = $seqin->next_seq) {
    $seqout->write_seq($seq);
}
