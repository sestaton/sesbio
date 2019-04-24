#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;
use autodie qw(open);
use Bio::DB::HTS::Kseq;

my $usage = "perl $0 genome.fasta\n";
my $infile = shift or die $usage;

my ($n_count, $total, $seq_ct, $n_ratio, $n_ratio_t, $n_perc, $n_perc_t, $gap_tot) = (0, 0, 0, 0, 0, 0, 0, 0);

my $kseq = Bio::DB::HTS::Kseq->new($infile);
my $iter = $kseq->iterator;

my $gaplen = 73;
my $gapstr = 'N' x $gaplen;

say "Seq_ID\tLength\tfake-gap_count\tfake-gap_ratio\tfake-gap_perc";

while (my $seqobj = $iter->next_seq) { 
    my $id = $seqobj->name;
    my $seq = $seqobj->seq;

    my $seqlength = length($seq);
    if ($seqlength > 0) {
        $seq_ct++;
        $total += $seqlength;
	$n_count = () = $seq =~ /[catg]$gapstr[catg]/ig;
	$gap_tot = $n_count * $gaplen; # 73 is the artificial gap size
        $n_perc  = sprintf("%.2f",$gap_tot/$seqlength);
        $n_ratio = sprintf("%.2f",$gap_tot/($seqlength - $gap_tot));
        $n_ratio_t += $n_ratio;
        $n_perc_t  += $n_perc; 
        say join "\t", $name, $seqlength, $n_count, $n_ratio, $n_perc;
    }
}
close $in;

my $len_ave     = sprintf("%.0f",$total/$seq_ct);
my $n_ratio_ave = sprintf("%.2f",$n_ratio_t/$seq_ct);
my $n_perc_ave  = sprintf("%.2f",$n_perc_t/$seq_ct);

say "\nAverage length, N-ratio and N-percent for $infile:   $len_ave\t$n_ratio_ave\t$n_perc_ave";
