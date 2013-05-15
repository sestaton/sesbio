#!/usr/bin/env perl

use v5.12;
use strict;
use warnings;
use autodie qw(open);
use Data::Dump qw(dd dump);
use List::Util qw(sum);
use Getopt::Long;

my $usage = "$0 -i annot.tsv -c cls -o out -n seqct\n";

my $infile;
my $cls;
my $seqct;
my $outfile;
my %annot;
my %sfam_hitct;
my %fam_readct;
my $total_ct = 0;

GetOptions(
           'i|infile=s'    => \$infile,
           'o|outfile=s'   => \$outfile,
           'n|seqnum=i'    => \$seqct,
           'c|cls=s'       => \$cls,
           );

die $usage if !$infile; # or !$outfile or !$cls;

#my $cls_tot = get_cls_tot($cls);
#my ($seqs, $seqct) = fas2hash($fas_file);
#my $repfrac = $cls_tot / $seqct;
open my $in, '<', $infile;

while (<$in>) {
    chomp;
    my ($th, $ct) = split;
    if ($th =~ /((^RL.\-\w+)\-\d+)/) {
	$th = $2;
	#say join "\t", $th, $1, $2;
	$annot{$th} += $ct;
    }
    else {
	$annot{$th} += $ct;
    }
}
my $sum = sum values %annot;
for my $k (reverse sort { $annot{$a} <=> $annot{$b} } keys %annot) {
    my $perc = sprintf("%.2f", $annot{$k} / $sum);
    say join "\t", $k, $annot{$k}, $perc;
}
