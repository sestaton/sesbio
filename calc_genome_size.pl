#!/usr/bin/env perl

use 5.014;
use strict;
use warnings;
use autodie qw(open);
use Bio::Kseq;
use List::MoreUtils qw(uniq);
use List::Utils qw(sum);
use Statistics::Descriptive;
use Getopt::Long;

my $usage = "perl $0 -i infile -q query -t target -ml mer_len -fl filter_len [-e] [-pid]\n";
my $infile;
my $query; 
my $target;
my $mer_len;
my $filter_len;
my $evalue;
my $pid;

GetOptions(
	   'i|infile=s'       => \$infile,
	   'q|query=s'        => \$query,
	   't|target=s'       => \$target,
	   'ml|merlen=i'      => \$mer_len,
	   'fl|filter_len=i'  => \$filter_len,
	   'e|evalue=f'       => \$evalue,
	   'pid|perc_ident=i' => \$pid,
   );

die $usage if !$infile or !$query or !$target or !$mer_len;

$filter_len //= 50;
$evalue //= 10;
$pid //= 80;

my ($query_store, $qseq_ct, $query_bases) = store_seq_len($query);
say STDERR "Finished storing $$qseq_ct sequences in $query.";
my ($target_store, $tseq_ct, $target_bases) = store_seq_len($target);
#my $target_bases = sum values %$target_store;
say STDERR "Finished storing $$tseq_ct sequences in $target ($target_bases bases).";

my %matches;

open my $in, '<', $infile;

while (<$in>) {
    chomp;
    my ($qid, $sid, $E, $N, $Sprime, $S, $alignlen, $nident, $npos, $nmism, $pcident, $pcpos, $qgaps, 
	$qgaplen, $sgaps, $sgaplen, $qframe, $qstart, $qend, $sframe, $sstart, $send) = split;
    # qid      | query sequence identifier
    # sid      | subject (database) sequence identifier
    # E        | the expectation or E-value
    # N        | the number of scores considered jointly in computing E
    # Sprime   | the normalized alignment score, expressed in units of bits
    # S        | the raw alignment score
    # alignlen | the overall length of the alignment including any gaps
    # nident   | the number of identical letter pairs
    # npos     | the number of letter pairs contributing a positive score
    # nmism    | the number of mismatched letter pairs
    # pcident  | percent identity over the alignment length (as a fraction of alignlen)
    # pcpos    | percent positive letter pairs over the alignment length (as a fraction of alignlen)
    # qgaps    | number of gaps in the query sequence
    # qgaplen  | total length of all gaps in the query sequence
    # sgaps    | number of gaps in the subject sequence
    # sgaplen  | total length of all gaps in the subject sequence
    # qframe   | the reading frame in the query sequence (+0 for protein sequences in BLASTP and TBLASTN searches)
    # qstart   | the starting coordinate of the alignment in the query sequence
    # qend     | the ending coordinate of the alignment in the query sequence
    # sframe   | the reading frame in the subject sequence (+0 for protein sequences in BLASTP and BLASTX searches)
    # sstart   | the starting coordinate of the alignment in the subject sequence
    # send     | the ending coordinate of the alignment in the subject sequence
    if ($alignlen >= $filter_len && $pcident >= $pid) {
	if (exists $matches{$qid}) {
	    push @{$matches{$qid}}, $sid;
	}
	else {
	    $matches{$qid} = [ $sid ];
	}
    }
}
say STDERR "Finished summarizing mapped hits to each transcript.";

my %counts;
for my $k (keys %matches) { # this will calculate the number of hits per transcript
    #$counts{$_}++ for @{$matches{$k}};
    my $count = uniq(@{$matches{$k}});
    $counts{$k} = $count;
}
%matches = ();
say STDERR "Finished calculating the number of mapped reads per transcript.";

my %trans_cov;
for my $transc (keys %counts) {
    if (exists $store->{$transc}) {
	my $len_mapped = $counts{$transc} * 59;
	my $trans_len = $store->{$transc};
	my $trans_cov = $len_mapped / $trans_len;
	#my $trans_cov_rnd = int($trans_cov + $trans_cov/abs($trans_cov*2)); # round num
	my $trans_cov_nornd = sprintf("%.2f", $trans_cov);
	$trans_cov{$transc} = $trans_cov_nornd;
    }
}
%counts = ();
say STDERR "Finished calculating the per transcript coverage. Now writing output.";

my $stat = Statistics::Descriptive::Full->new();
for my $trans (reverse sort { $trans_cov{$a} <=> $trans_cov{$b} } keys %trans_cov) {
    #say join "\t", $trans, $trans_cov{$trans};
    $stat->add_data($trans_cov{$trans});
}

my $count = $stat->count;
my $mean  = $stat->mean;
my $medi  = $stat->median;
my $var   = $stat->variance;
my $min   = $stat->min;
my $max   = $stat->max;
my $sd    = $stat->standard_deviation;

say STDERR join "\t", "Count", "Mean", "Median", "Variance", "SD", "Min", "Max";
say STDERR join "\t", $count, $mean, $medi, $var, $sd, $min, $max;
say STDERR join "\t", "cval_by_mean", "cval_by_median", "cval_by_min", "cval_by_max";
say STDERR join "\t", ($target_bases/$mean), ($target_bases/$medi), ($target_bases/$min), ($target_bases/$max);

undef $stat;
exit;
#
# subs
#
sub store_seq_len {
    my ($query) = @_;

    my %store;
    my $kseq = Bio::Kseq->new($query);
    my $it = $kseq->iterator;

    my ($n, $seq_ct, $base_ct) = (0, 0, 0);
    while( my $seq = $it->next_seq) {
	$seq_ct++;
	$base_ct += length($seq->{seq});
	$store{ $seq->{name} } = length($seq->{seq});
    }
    return \%store, \$seq_ct, \$base_ct;
}
