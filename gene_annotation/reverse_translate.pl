#!/usr/bin/env perl

use strict;
use warnings;
use Bio::Tools::SeqPattern;
use Bio::SeqIO;
use Bio::PrimarySeq;
use Data::Dumper;

my $usage = "\n$0 pepfile ntfile\n\n";
my $pep = shift or die $usage;
my $nt = shift or die $usage;

open my $out, '>', $nt or die "\nERROR: Could not open file: $nt\n";

my $myCodonTable = Bio::Tools::CodonTable->new();

my $seq_in = Bio::SeqIO->new(-file => $pep,
			     -format => 'fasta');

while( my $seq = $seq_in->next_seq ) {
    my $seq_id = $seq->id;
    my $seq_pep = $seq->seq;
   
    my $pattern = Bio::Tools::SeqPattern->new(-SEQ =>$seq_pep, -TYPE =>'Amino');
    my $nuc = $pattern->backtranslate;
    print $out ">$seq_id\n$nuc->{str}\n";
}

close $out;
