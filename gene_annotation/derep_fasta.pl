#!/usr/bin/perl -w

# remove duplicate entries from a fasta file
# 4/12/12 SES

use strict;
use Bio::SeqIO;
use Getopt::Long;

my $usage = "\n$0 -i in_fas -o derep_fas\n\n";
my $infile;
my $outfile;

GetOptions(
	   'i|infile=s'  => \$infile,
	   'o|outfile=s' => \$outfile,
	   );

die $usage if !$infile || !$outfile;

my $seq_in  = Bio::SeqIO->new(-file => $infile, -format => 'fasta');
my $seq_out = Bio::SeqIO->new(-file => ">$outfile", -format => 'fasta');

my %seqhash;

while(my $seq = $seq_in->next_seq) {
    unless(exists $seqhash{$seq->id}) {
	$seq_out->write_seq($seq);
    }
    $seqhash{$seq->id} = $seq->seq;
}

exit;
