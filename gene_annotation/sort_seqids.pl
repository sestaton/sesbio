#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;
use Bio::SeqIO;
use Getopt::Long;
use Time::HiRes qw(gettimeofday);

my $infile;
my $outfile;
my $usage = "$0 -i in -o sorted\n";

GetOptions(
	   'i|infile=s'   => \$infile,
	   'o|outfile=s'  => \$outfile,
	   );

if (!$infile or !$outfile) {
    say $usage,"\nERROR: No input was given. Exiting.";
    exit(1);
}

my $seq_in = Bio::SeqIO->new(-file => $infile, -format => 'fasta');
open my $out, ">", $outfile or die "\nERROR: Could not open file: $outfile";

my %seqhash;
my $seqct = 0;
my $t0 = gettimeofday();

while (my $seq = $seq_in->next_seq()) {
    $seqhash{$seq->id} = $seq->seq;
    $seqct++;
}

for my $k (sort keys %seqhash) {
    say $out join "\n", ">".$k,$seqhash{$k};
}
close $out;

my $t1 = gettimeofday();
my $elapsed = $t1 - $t0;
my $time = sprintf("%.2f",$elapsed);

say "========== Sorted $seqct sequences in $time seconds.";
