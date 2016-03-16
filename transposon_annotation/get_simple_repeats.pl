#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;
use Bio::SeqIO;

my $usage = "$0 seq.fas";
my $fasta = shift or die $usage;

#my $seqio = Bio::SeqIO->new( -fasta => 'fasta', -file => $fasta );

my %di = ('AA' => 0, 'AC' => 0, 
	  'AG' => 0, 'AT' => 0, 
	  'CA' => 0, 'CC' => 0, 
	  'CG' => 0, 'CT' => 0, 
	  'GA' => 0, 'GC' => 0, 
	  'GG' => 0, 'GT' => 0, 
	  'TA' => 0, 'TC' => 0, 
	  'TG' => 0, 'TT' => 0);
my %mono = ('A' => 0, 'C' => 0, 'G' => 0, 'T' => 0);

#my $seqct  = 0;
#my $repct  = 0;
#my $dict   = 0;
#my $monoct = 0;
#my $diratio;
#my $monoratio;

my $hastot = 0;
for my $repeat_ratio (qw(0.30 0.40 0.50 0.60 0.70 0.80 0.90)) {
    my $seqct  = 0;
    my $repct  = 0;
    my $dict   = 0;
    my $monoct = 0;
    my $diratio;
    my $monoratio;

    say "working on $repeat_ratio";
    my $seqio = Bio::SeqIO->new( -fasta => 'fasta', -file => $fasta );

    while (my $seqin = $seqio->next_seq) {
	$seqct++;
	my $seq = $seqin->seq;
	my $len = $seqin->length;
	
	my $counted = 0;
	for my $mononuc (keys %mono) {
	    while ($seq =~ /$mononuc/ig) { $monoct++ };
	    $monoratio = sprintf("%.2f", $monoct/$len);
	    if ($monoratio >= $repeat_ratio) {
		$counted++;
		$repct++ unless $counted > 1;	
	    }
	    $monoct = 0;
	}

	for my $dinuc (keys %di) {
	    while ($seq =~ /$dinuc/ig) { $dict++ };
	    $diratio = sprintf("%.2f", $dict*2/$len);
	    if ($diratio >= $repeat_ratio) {
		$counted++;
		$repct++ unless $counted > 1;
	    }
	    $dict = 0;
	}
	#$counted = 0;
    }

    my $repfrac = sprintf("%.2f", $repct/$seqct*100);
    say STDERR "Total sequences                             : $seqct" unless $hastot;
    say STDERR "Num simple repeats filtered at $repeat_ratio: $repct";
    say STDERR "Repeat fraction at $repeat_ratio            : $repfrac";
    $hastot  = 1;
    $repct   = 0;
    $repfrac = 0;
    $seqct   = 0;
}
