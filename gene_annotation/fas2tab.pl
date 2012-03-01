#!/usr/bin/perl -w

use strict;
use Bio::SeqIO;

my $usage = "$0 infile > out\n";
my $infile = shift or die $usage;

my $seq_in = Bio::SeqIO->new(-file => $infile,
			     -format => 'fasta');

#open( my $in, '<', $infile) or die "\nERROR: Could not open file: $infile\n";

my $header;
my @seq;

#while(my $line = <$in>) {
while(my $seq = $seq_in->next_seq) {
    #chomp $line;
    #if ($line =~ m/^\>/) {
	#$header = $line;
    #} else {
	#push(@seq, $line);
    #}
    #my $fa = join('',@seq);
    #$fa =~ s/\s//g;
    print join("\t",($seq->id,$seq->seq)),"\n";
    #print "$header\t$fa\n";
}

#close($in);

exit;
