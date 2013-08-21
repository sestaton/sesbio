#!/usr/bin/env perl

# NB: The input files need to be sorted -- that is, each pair needs to be in the
#     same order. Also, this is probably too slow to be useful.

use strict;
use warnings;
use Bio::SeqIO;
use Getopt::Long;

my $usage = "\nUSAGE: $0 -p1 pair1.fastq -p2 pair2.fastq -o outfile\n\n";
my $pair1;
my $pair2;
my $outfile;
my $variant;

GetOptions(
           'p1|pair1=s'    => \$pair1,
	   'p2|pair2=s'    => \$pair2,
	   'o|outfile=s'   => \$outfile,
           'v|variant=s'   => \$variant,
          );

# open the infile or die with a usage statement
if ($pair1 && $pair2 && $outfile)  {
    print "\n=========== Creating interleaved fastq file $outfile ...\n";
} 
else {
    die $usage;
}

$variant = defined($variant) ? $variant : 'sanger';

# create SeqIO objects to read in and to write outfiles
my %seq_in = ( 
               'pair_1' => Bio::SeqIO->new('-file'    => "<$pair1",
					   '-variant' => $variant,
                                           '-format'  => 'fastq'),
               'pair_2' => Bio::SeqIO->new('-file'    => "<$pair2",
					   '-variant' => $variant,
					   '-format'  => 'fastq'),
             );

my $seq_out = Bio::SeqIO->new(-format   => 'fastq', 
			      -variant  => $variant, # this would need to be changed for Illumina Fastq files
			      -file     => ">$outfile");

while (my $seq1 = $seq_in{'pair_1'}->next_seq) {
    while (my $seq2 = $seq_in{'pair_2'}->next_seq) {
	if ($seq1->id =~ /\/1$|\ :1/ && $seq2->id =~ /\/2$|\ :2/) {               
	    $seq_out->write_seq($seq1,$seq2);
	}
    }
}
    
    
