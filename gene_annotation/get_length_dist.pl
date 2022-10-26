#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;
use File::Basename;
use Bio::DB::HTS::Kseq;
use Getopt::Long;

my $infile;
my $outfile;

GetOptions(# Required arguments
           'i|infile=s'         => \$infile,
           'o|outfile=s'        => \$outfile,
           );

# open the infile or die with a usage statement
usage() and exit(1) if !$infile or !$outfile;
open my $out, '>', $outfile or die "\nERROR: Could not open file: $outfile\n";

if (! -e $infile) {
    say STDERR "\n[ERROR]: '$infile' does not exist. Exiting.\n";
    exit(1);
}

my $kseq = Bio::DB::HTS::Kseq->new($infile);
my $iter = $kseq->iterator();

my ($nreads, $bases) = (0, 0, 0);
my ($min, $max);
my %dist;

while (my $seqobj = $iter->next_seq) {
    $nreads++;
    my $id = $seqobj->name;
    my $seq = $seqobj->seq;
    my $len = length($seq);

    if ($len < 1) {
        say STDERR "\n[ERROR]: '$id': len: $len\n";
        exit(1);
    }

    $bases += $len;
    $min = defined $min && $min < $len ? $min : $len;
    $max = defined $max && $max > $len ? $max : $len;
    $dist{$len}++;
}

my $mean_len = sprintf("%.0f", $bases / $nreads);
say STDERR join "\t", "Filename", "Min", "Max", "Mean";
say STDERR join "\t", $infile, $min, $max, $mean_len;

say $out join "\t", "Length", "Count";
for my $len (sort { $a <=> $b } keys %dist) {
    say $out join "\t", $len, $dist{$len};
}
close $out;

exit;
#
# sub
#
sub usage {
    my $script = basename($0);
    print STDERR <<END
USAGE: $0 -i seqs.fasta -o outreport

This script takes as input a FASTA/Q and outputs a length distribution such as:

Length  Count
1009      1
1038      1
1134      1
1152      1

NB: The lengths will report in sort order (numeric, increasing).

END
}
