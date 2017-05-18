#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;
use autodie;
use Sort::Naturally;
use Bio::DB::HTS::Kseq;
use Bio::DB::HTS::Faidx;
#use Data::Dump::Color;

my $usage = "$0 seq";
my $seqfile = shift or die $usage;

my $windows = make_windows($seqfile);
#dd $windows and exit(1);

my $index = Bio::DB::HTS::Faidx->new($seqfile);

my %bins;
my $binct = 0;
for my $chr (nsort keys %$windows) {
    for my $window (nsort @{$windows->{$chr}}) {
	my ($wct, $start, $end) = split /\|\|/, $window;

	$start = $start == 0 ? $start+1 : $start;
	my $id = "$chr:$start-$end";
	my ($seq, $len)  = $index->get_sequence($id);
	unless ($len) {
	    say STDERR "DEBUG:\n$id has zero length.";
	    exit;
	}
	my $gc_num = $seq =~ tr/GCgc/GCgc/;
	my $gc_content = sprintf("%.2f", $gc_num/$len);
	say join "\t", $chr, $start, $end, $gc_content;
    }
}

sub make_windows {
    my ($sequence) = @_;

    my $wsize = 1e6;

    my %windows;
    my $kseq = Bio::DB::HTS::Kseq->new($sequence);
    my $iter = $kseq->iterator;

    while (my $seqobj = $iter->next_seq) {
	my $seq = $seqobj->{seq};
	my $id  = $seqobj->{name};
	my $length = length($seq);
	my $remainder = $length;
        my ($total, $start, $end, $wct, $chunk_size) = (0, 0, 0, 0, 0);
        my $steps = sprintf("%.0f", $length/$wsize);

	$steps = $length % $wsize == 0 ? $steps : $steps + 1;
        $steps = 0 if $length <= $wsize;

	for (0..$steps) {
	    last if $remainder == 0;
            if ($remainder < $wsize) {
		$start = $end;
		$end = $length;
		$remainder = 0;
	    }
	    else {
		if ($wct) {
		    $start = $end;
		    $end = $start + $wsize;
		    $remainder -= $wsize;
		}
		else {
		    $start = 0;
		    $end = $wsize;
		    $remainder -= $wsize;
		}
	    }
	    push @{$windows{ $seqobj->{name} }}, join "||", $wct, $start, $end;
	    $wct++;
	}
    }
    
    return \%windows;
}

