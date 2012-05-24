#!/usr/bin/perl -w

use strict;
use Bio::AlignIO;

my $usage = "USAGE: cnv_alignformat.pl infile outfile\n"; 
my $inputfilename = $ARGV[0] or die $usage;
my $outputfilename = $ARGV[1] or die $usage;

my $alignin  = Bio::AlignIO->new(-file => $inputfilename , '-format' => 'fasta');
my $alignout = Bio::AlignIO->new(-file => ">$outputfilename" , '-format' => 'clustalw');
    # note: we quote -format to keep older perl's from complaining.

while ( my $aln = $alignin->next_aln() ) {
    $alignout->write_aln($aln);
}
