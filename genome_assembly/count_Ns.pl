#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;
use autodie qw(open);
use Bio::DB::HTS::Kseq;

my $usage = "perl $0 seqfile.fastq\n";
my $infile = shift or die $usage;

my $knseq = Bio::DB::HTS::Kseq->new($infile);
my $nt_it = $knseq->iterator;

my ($n_count, $total, $seq_ct, $n_ratio, $n_ratio_t, $n_perc, $n_perc_t) = (0, 0, 0, 0, 0, 0, 0);

say "Seq_ID\tLength\tN_count\tN_ratio\tN_perc";

while (my $seq_in = $nt_it->next_seq) {
    my $name = $seq_in->name;
    my $seq = $seq_in->seq;
    my $seqlength = length($seq);
    if ($seqlength > 0) {
        $seq_ct++;
        $total += $seqlength;
        $n_count = ($seq =~ tr/Nn//);
        $n_perc  = sprintf("%.2f",($n_count/$seqlength)*100);
        $n_ratio = sprintf("%.2f",$n_count/($seqlength - $n_count));
        $n_ratio_t += $n_ratio;
        $n_perc_t  += $n_perc; 
        say join "\t", $name, $seqlength, $n_count, $n_ratio, $n_perc;
    }
}

my $len_ave     = sprintf("%.0f",$total/$seq_ct);
my $n_ratio_ave = sprintf("%.2f",$n_ratio_t/$seq_ct);
my $n_perc_ave  = sprintf("%.2f",$n_perc_t/$seq_ct);

say "\nAverage length, N-ratio and N-percent for $infile:   $len_ave\t$n_ratio_ave\t$n_perc_ave";
