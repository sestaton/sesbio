#!/usr/bin/env perl

## written to transfer annotations to new assemblies using RATT from EMBL

use strict;
use warnings;
use Bio::SeqIO;
use Getopt::Long;

my $usage = "perl $0 -i seqs_anno.gb -o seqs_anno.embl\n";
my $infile;
my $outfile;

GetOptions(
           'i|infile=s'  => \$infile,
           'o|outfile=s' => \$outfile,
           );

die $usage if !$infile || !$outfile;

my $seqin  = Bio::SeqIO->new(-file => $infile,     -format => 'genbank');
my $seqout = Bio::SeqIO->new(-file => ">$outfile", -format => 'embl');

while (my $seq = $seqin->next_seq) {
    $seqout->write_seq($seq);
}
