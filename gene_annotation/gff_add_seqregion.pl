#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;
use Sort::Naturally;
use Bio::DB::HTS::Kseq;
use Getopt::Long;

my %opts;
GetOptions(\%opts, 'genome|g=s', 'gff|f=s', 'outfile|o=s');
if (!%opts) {
    say "\nERROR: Command line not parsed correctly.\n";
    say "USAGE: perl add_seqregion.pl -g genome.fas -f genome.gff3 -o genome_seqreg-added.gff3\n";
    exit 1;
}

my %lens;
my $kseq = Bio::DB::HTS::Kseq->new($opts{genome});
my $iter = $kseq->iterator;
while (my $seqobj = $iter->next_seq) {
    my $id = $seqobj->name;
    my $seq = $seqobj->seq;
    my $len = length($seq);
    $lens{$id} = $len;
}

open my $out, '>', $opts{outfile} or die "\nERROR: Could not open file: $opts{outfile}\n";
say $out '##gff-version 3';

for my $seq (nsort keys %lens) {
    say $out join q{ }, '##sequence-region', $seq, '1', $lens{$seq};
}

open my $in, '<', $opts{gff} or die "\nERROR: Could not open file: $opts{gff}\n";

while (my $line = <$in>) {
    chomp $line;
    next unless $line =~ /\S/;
    next if $line =~ /^#/;
    say $out $line;
}

close $in;
close $out;
