#!/usr/bin/perl -w

# get_seq_length.pl


use strict;
use Bio::SeqIO;
       
#----------------+
# VARIABLE SCOPE |
#----------------+
# get command-line arguments
my $usage = "\nUSAGE: $0 infile outfile\n\nThe outfile will contain a summary of all the contigs.\n\n";
my $infile = $ARGV[0] || die "\n","ERROR: No infile was given at the command line\n",$usage;
my $outfile = $ARGV[1] || die "\n","ERROR: No outfile was given at the command line\n",$usage;

# create SeqIO objects to read in and to write outfiles
my $seq_in  = Bio::SeqIO->new( -format => 'fasta', 
                               -file => $infile); 

open( OUT, ">$outfile") || die "\nERROR: Could not open file: $outfile\n";
while(my $seq = $seq_in->next_seq()) {
    print OUT $seq->id."\t".$seq->length."\n";
}

close(OUT);
exit;
