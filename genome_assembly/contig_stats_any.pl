#!/usr/bin/env perl

# chienchi@lanl.gov
# generate some stats for Newbler or Velvet (>0.7.6) contigs.
# Assume contigs fasta file is in the assembly output folder
#   which contains other files, like Log (velvet), stats.txt(velvet),
#   and 454ContigGraph.txt(Newbler)
# 20100423
# no flag -t can generate some statistics with contig file only.
# 20100708
#
# modified by S. Evan Staton <statonse at gmail dot com>

use 5.010;
use strict;
use warnings;
use File::Basename;
use Getopt::Long;

GetOptions("help|?"   => sub {Usage()});

if (scalar(@ARGV) != 1){&Usage;}

my ($file_name, $path, $suffix)=fileparse("$ARGV[0]", qr/\.[^.]*/);
my ($len,$total)=(0,0);
my @x;
my $seq_num;
my $seq;
my $reads_num;
my $kmer_cov;
my $GC_num;
my $GC_content;
my $id_to_reads;
my $id_to_cov;
my $used_percent;
my $singleton;
my $exp_cov;
my ($over100k_bases,$over50k_bases,$over25k_bases,$over10k_bases,$over5k_bases,$over3k_bases,$over2k_bases,$over1k_bases)=(0,0,0,0,0,0,0,0,0);
my ($over100k_reads,$over50k_reads,$over25k_reads,$over10k_reads,$over5k_reads,$over3k_reads,$over2k_reads,$over1k_reads)=(0,0,0,0,0,0,0,0,0);

while (<>) {
    chomp;
    if (/^[\>\@]/) {
	$seq_num++;
	if ($len > 0) {
	    stats($len);
	    push @x, $len;
	    $GC_num = $seq =~ tr/GCgc/GCgc/;
	    $GC_content = $GC_num/$len;
	}
	$len = 0;
	$seq = "";
    }
    else {
	s/\s//g;
	$len += length($_);
	$seq .= $_;
    }
}

if ($len > 0) {
    stats($len);
    push @x, $len;
    $GC_num = $seq =~ tr/GCgc/GCgc/;
    $GC_content = $GC_num/$len;
}
@x= sort { $b <=> $a } @x;
my $N50;
my $N90;
my ($count,$half)=(0,0);
for (my $j = 0; $j < @x; $j++) {
    $count+=$x[$j];
    if (($count >= $total / 2) && ($half == 0)){
	$N50=$x[$j];
	$half=$x[$j]
    } elsif ($count >= $total * 0.9){
	$N90=$x[$j];
	last;
    }
}

my ($top10, $top20, $top40, $top100);
for (0..99) {
    $top10 += $x[$_] if ($_ < 9 and $x[$_]);
    $top20 += $x[$_] if ($_ < 19 and $x[$_]);
    $top40 += $x[$_] if ($_ < 39 and $x[$_]);
    $top100 += $x[$_] if ($_ < 99 and $x[$_]);
}

say "Contigs_number:\t$seq_num";
say "N50:\t$N50";
say "N90:\t$N90";
say "Max:\t$x[0]";
say "Min:\t$x[-1]";
say "Total_bases:\t$total";
say "Top10_bases:\t$top10";
say "Top20_bases:\t$top20";
say "Top40_bases:\t$top40";
say "Top100_bases:\t$top100";
say ">100kb_bases:\t$over100k_bases";
say ">50kb_bases:\t$over50k_bases";
say ">25kb_bases:\t$over25k_bases";
say ">10kb_bases:\t$over10k_bases";
say ">5kb_bases:\t$over5k_bases";
say ">3kb_bases:\t$over3k_bases";
say ">2kb_bases:\t$over2k_bases";
say ">1kb_bases:\t$over1k_bases";

#
# Subs
#
sub Usage {
    say STDERR "perl $0 <contigs.fasta>";
    say STDERR "     -help       print usage";
    exit;
}

sub stats {
    my $len = shift;
    $total += $len;
    $over100k_bases += $len if ($len > 100000);
    $over50k_bases += $len if ($len > 50000);
    $over25k_bases += $len if ($len > 25000);
    $over10k_bases += $len if ($len > 10000);
    $over5k_bases += $len if ($len > 5000);
    $over3k_bases += $len if ($len > 3000);
    $over2k_bases += $len if ($len > 2000);
    $over1k_bases += $len if ($len > 1000);
}
