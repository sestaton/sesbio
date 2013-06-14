#!/usr/bin/perl -w

use strict;
#use lib qw(/iob_home/jmblab/statonse/apps/perlmod/bioperl-run/lib);
#use Bio::Tools::CodonTable;
use Bio::Tools::SeqPattern;
use Bio::SeqIO;
use Bio::PrimarySeq;
use Data::Dumper;

my $usage = "\n$0 pepfile ntfile\n\n";
my $pep = shift or die $usage;
my $nt = shift or die $usage;

open(my $out, '>', $nt) or die "\nERROR: Could not open file: $nt\n";

my $myCodonTable   = Bio::Tools::CodonTable->new();

my $seq_in = Bio::SeqIO->new(-file => $pep,
			     -format => 'fasta');

#my $seq_out = Bio::SeqIO->new(-file => ">$nt",
#			      -format => 'fasta');

while( my $seq = $seq_in->next_seq ) {
    my $seq_id = $seq->id;
    my $seq_pep = $seq->seq;
   
    my $pattern = Bio::Tools::SeqPattern->new(-SEQ =>$seq_pep, -TYPE =>'Amino');
    #my $prim_seq = Bio::PrimarySeq->new(-seq => $seq_pep);
    #my $iupac_str = $myCodonTable->reverse_translate_all($prim_seq);
    my $nuc = $pattern->backtranslate;
    #my $prim_seq = Bio::PrimarySeq->new(-seq => $nuc);
    #$seq_out->write_seq($prim_seq);
    #print Dumper $nuc;
    print $out ">$seq_id\n$nuc->{str}\n";
}

#close($out);

exit;
