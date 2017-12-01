#!/usr/bin/env perl

## ABOUT: This is a script that takes a TBLASTN report (-outfmt 6) of protein alignments to a set of full-length TEs,
##        and will remove overlapping, redundant alignments and also report the strand of the query and subject alignments (very
##        useful to interpreting the potential role of the acquired gene products).
##
## NB: As of now, there are hard-coded thresholds for 50% identity and 50% coverage of the query alignment to the subject (TE).
##     This may be changed below but cautious if you reduce these thresholds as it will increase the error rate.
##
##     Last, it is imperative to filter non-specific proteins, binding domains, and TE-related proteins from the report. 
##     Otherwise, the amount of gene fragments acquired will be inflated. I may add this function in the future (currently handled
##     by a separate script).
##
## AUTHOR: S. Evan Staton <evan at evanstaton.com>
## LINK: https://github/sestaton/sesbio
## LICENSE: MIT @ S. Evan Staton, 2017

use 5.010;
use strict;
use warnings;
use autodie;
use Sort::Naturally;
use Bio::DB::HTS::Kseq;
use Data::Dump::Color;

my $usage = "query_protfile.faa blast.bln";
my $fasta = shift or die $usage;
my $blast = shift or die $usage;

my $lens = get_lengths($fasta);

open my $in, '<', $blast;

say join "\t", "qid", "sid", "pid", "alnlen", "mism", "gaps", "qstart", "qend", "sstart", "send", "evalue", "bits", "qstrand", "sstrand";
my %hits;
while (<$in>) {
    chomp;
    my @f = split /\t/;
    
    # adjust regex, which will inevitably be unique for nearly every genome
    ## TODO: print warning about ID formats, which need to be handling carefully for this to work
    my ($chr, $sstart, $send) = ($f[1] =~ /(MtrunA17(?:MT|CP)?(?:Chr\d+)?(?:\w+\d+)?)_(\d+)_(\d+)(?:_fragment_)?/);
    unless (defined $sstart && defined $send) {
	say "\nERROR: $_";
	exit;
    }

    my $slen = $send - $sstart + 1;
    my $qlen = $lens->{$f[0]};
    defined $qlen && $qlen ne '' 
	or die "\nERROR: could not get length for $f[0]";

    my ($qstrand, $sstrand, $qaln_len, $saln_len, $qcoords, $scoords);
    my ($qaln_start, $qaln_end,$saln_start, $saln_end);
    if ($f[8] > $f[9]) {
	$saln_len = $f[8] - $f[9] + 1;
	$sstrand = '-';
	($saln_start, $saln_end) = @f[9,8];
	$scoords = join "||", @f[9,8];
    }
    else { 
	$saln_len = $f[9] - $f[8] + 1;
	$sstrand = '+';
	($saln_start, $saln_end) = @f[8,9];
	$scoords = join"||", @f[8,9];
    }
    if ($f[6] > $f[7]) {
	$qaln_len = $f[6] - $f[7] + 1;
	$qstrand = '-';
	($qaln_start, $qaln_end) = @f[7,6];
	$qcoords = join"||", @f[7,6];
    }
    else {
	$qaln_len = $f[7] - $f[6] + 1;
	$qstrand = '+';
	($qaln_start, $qaln_end) = @f[6,7];
	$qcoords = join"||", @f[6,7];
    }
 
    my $qaln_perc = sprintf("%.2f",($qaln_len/$qlen)*100);
    my $saln_perc = sprintf("%.2f",($saln_len/$slen)*100);

    if ($qaln_perc >= 50 && $f[2] >= 50) { 
	push @{$hits{$f[1]}}, join "~~", @f[0..5], $qaln_start, $qaln_end, $saln_start, $saln_end, @f[10..11], $qstrand, $sstrand;
    }
}

# this works by comparing alignments to each TE, sorted by coordinate
for my $hit (nsort keys %hits) {
    # get this first element by coordinate
    my ($first) = (map { $_->[0] } sort { $a->[9] <=> $b->[9] } map { [ $_, split /\~\~/ ] } @{$hits{$hit}})[0];

    say join "\t", split /\~\~/, $first;

    # coordinates of alignments to first element, sorted by coordinate
    my ($prev_start, $prev_end) = (split /\~\~/, $first)[8,9];
    for my $aln (@{$hits{$hit}}) {
	next if $aln eq $first;
	my ($this_start, $this_end) = (split /\~\~/, $aln)[8,9];
	if ($this_start > $prev_end) { # and $this_end > $prev_end) {
	    say join "\t", split /\~\~/, $aln;
	    ($prev_start, $prev_end) = ($this_start, $this_end);
	}
    }
}

exit;
#
# methods
#
sub get_lengths {
    my ($fas) = @_;

    my %lens;
    my $kseq = Bio::DB::HTS::Kseq->new($fas);
    my $iter = $kseq->iterator;

    while (my $obj = $iter->next_seq) {
        my $id = $obj->name;
        my $len =length($obj->seq);

        if ($len > 0) { 
            $lens{ $id } = $len;
        }
        else {
            die "ERROR: zero length sequence for: $id\n";
        }
    }

    return \%lens;
}
