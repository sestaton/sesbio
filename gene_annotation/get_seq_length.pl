#!/usr/bin/env perl

# TODO: Add minimal POD
use 5.010;
use strict;
use warnings;
use Bio::SeqIO;
use Getopt::Long;
       
#
# Vars
#
my $usage = "\nUSAGE: $0 -i infile -o outfile\n\nThe outfile will contain a summary of all the contig lengths formatted as:\nread_name\tlength(bp)\n\n";
my $infile; 
my $outfile;

GetOptions(
          'i|infile=s'  => \$infile,
          'o|outfile=s' => \$outfile,
          );

if (!$infile || !$outfile) {
    die "\nERROR: Command line not parsed correctly. Exiting.\n",$usage;
}

# create SeqIO objects to read in and to write outfiles
my $seq_in  = Bio::SeqIO->new( -format => 'fasta', -file => $infile); 

open my $out, '>', $outfile or die "\nERROR: Could not open file: $outfile\n";
while(my $seq = $seq_in->next_seq()) {
    say join "\t", $out $seq->id, $seq->length;
}

close $out;
