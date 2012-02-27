#!/usr/bin/perl -w

use strict;
use lib qw(/usr/local/bioperl/latest/); 
use Bio::SeqIO;
use Getopt::Long;

my $usage = "USAGE: pairs2interleaved.pl -p1 pair1.fastq -p2 pair2.fastq -o outfile\n";
my $pair1;
my $pair2;
my $outfile;

GetOptions(
	    "p1|pair1=s"    => \$pair1,
	    "p2|pair2=s"    => \$pair2,
	    "o|outfile=s"   => \$outfile,
          );

# open the infile or die with a usage statement
if ($pair1 && $pair2 && $outfile )  {
    #open (OUTFILE, ">$outfile")
    print "\n=========== Creating interleaved fastq file $outfile ...\n";
} 
else {
    die "\n",$usage,"\n";
}

# create SeqIO objects to read in and to write outfiles
my %seq_in = ( 
               'pair_1' => Bio::SeqIO->new('-file'    => "<$pair1",
					   '-variant' => 'illumina',
                                           '-format'  => 'fastq'),
               'pair_2' => Bio::SeqIO->new('-file'    => "<$pair2",
					   '-variant' => 'illumina',
					   '-format'  => 'fastq'),
             );

my $seq_out = Bio::SeqIO->new(-format   => 'fastq', 
			      -variant  => 'illumina',
			      -file     => ">$outfile");

while (my $seq1 = $seq_in{'pair_1'}->next_seq) { #     && my $seq2 = $seq_in{'pair_2'}->next_seq) {

   
    #$seq_out->write_seq($seq1);
    
    while (my $seq2 = $seq_in{'pair_2'}->next_seq) {

	#next if ($seq2->seq =~ /^n+$/);
       
	#die "\nThese reads are out of order!\n" unless $seq1->id eq$seq2->id;
	#foreach ($seq1 && $seq2) {
	#    $seq_out->write_seq($seq1);
	
	$seq_out->write_seq($seq1,$seq2);
	   
	#}

	#die "\nThese reads are out of order!\n" unless $seq1->id eq $seq2->id;
    }
}
    
    
