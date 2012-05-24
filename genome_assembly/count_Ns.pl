#!/usr/bin/perl -w

use strict;
use Bio::SeqIO;

my $usage = "$0 infile\n";
my $in = shift or die $usage;

my $seq_in = Bio::SeqIO->new(-file => $in, -format => 'fasta');

my ($n_count, $total, $seq_ct, $n_ratio, $n_ratio_t, $n_perc, $n_perc_t) = (0, 0, 0, 0, 0, 0, 0);

print "Seq_ID\tLength\tN_count\tN_ratio\tN_perc\n";

while(my $seq = $seq_in->next_seq) {
    if ($seq->length > 0) {
        $seq_ct++;
        $total += $seq->length;
        $n_count = ($seq->seq =~ tr/Nn//);
        $n_perc  = sprintf("%.2f",$n_count/$seq->length);
        $n_ratio = sprintf("%.2f",$n_count/($seq->length - $n_count));
        $n_ratio_t += $n_ratio;
        $n_perc_t  += $n_perc; 
        print join("\t",($seq->id,$seq->length,$n_count,$n_ratio,$n_perc)),"\n";
    }
}

my $len_ave     = sprintf("%.0f",$total/$seq_ct);
my $n_ratio_ave = sprintf("%.2f",$n_ratio_t/$seq_ct);
my $n_perc_ave  = sprintf("%.2f",$n_perc_t/$seq_ct);

print "\nAverage length, N-ratio and N-percent for $in :   $len_ave\t$n_ratio_ave\t$n_perc_ave\n";
