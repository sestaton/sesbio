#!/usr/bin/perl -w
#_________________________________________________________________+
#                                                                 |
# cnv_genbank2fasta.pl - convert genbank to fasta                 |
#_________________________________________________________________+
#                                                                 |
# Description: Take a genbank file as input and write a fasta     |
# for the genbank record                                          |
#                                                                 |
# Author: S. Evan Staton                                          |
# Contact: statonse<at>uga.edu                                    |
# Started: 4.16.10                                                |
# Updated:                                                        |
#                                                                 |
# Usage: cnv_genbank2fasta.pl infile outfile                      |
#_________________________________________________________________+
# TODO: Make a separate fasta for each genbank entry (if it does 
# not already do this) - think about a way to write files for each 
# gene region or specify a region at the command-line (this is explained
# in Jason Stajich's bioperl tutorial)
#-----------+
# INCLUDES  |
#-----------+
use strict; 
use Bio::SeqIO;
       
#----------------+
# VARIABLE SCOPE |
#----------------+
# get command-line arguments, or die with a usage statement
my $usage = "cnv_genbank2fasta.pl infile outfile\n";
my $infile = shift or die $usage;
my $outfile = shift or die $usage;

# create SeqIO objects to read in and to write outfiles
my $seq_in = Bio::SeqIO->new(-format=>'genbank', -file=>"$infile");
my $seq_out = Bio::SeqIO->new(-format=>'fasta', -file=>">$outfile");

#Do it
my $seq;
while (my $seq = $seq_in->next_seq) {
  $seq_out->write_seq($seq);
}
exit;
