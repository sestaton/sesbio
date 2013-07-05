#!/usr/bin/env perl

use v5.10;
use strict;
use warnings;
use autodie qw(open);

my $usage = "perl $0 blast fas\n";
my $blast = shift or die $usage;
my $fasta = shift or die $usage;

open my $in, '<', $blast;
open my $fas, '<', $fasta;

my %match_range;

while (my $l = <$in>) {
    chomp $l;
    my @f = split "\t", $l;
    if (@f) { # check for blank lines in input
	next if exists $match_range{$f[0]};
	$match_range{$f[0]} = join "|", $f[6], $f[7];
    }
    else {
	#say "Here is the error at line: $.", join "\t", @f;
    }
}
close $in;

my ($scrSeqCt, $validscrSeqCt) = (0, 0);

{
    local $/ = '>';
    
    while (my $line = <$fas>) {
	$line =~ s/>//g;
	next if !length($line);
	my ($seqid, @seqs) = split /\n/, $line; 
	my $seq = join '', @seqs;
	my $useq = uc($seq);
	$scrSeqCt++ if defined $seq;
	$seqid =~ s/\s.*//;
	if (exists $match_range{$seqid}) {
	    my ($match_start, $match_end) = split /\|/, $match_range{$seqid};
	    if (defined $match_start && defined $match_end) {
		my $match_length = $match_end - $match_start;
		#say join "\t", $match_start, $match_end, $match_length;
		if ($match_length >= 50) {
		    $validscrSeqCt++;
		    my $seq_match = substr $useq, $match_start, $match_length;
		    say join "\n", ">".$seqid, $seq_match;
		}
	    }
	    else {
		#say "Not defined start: $match_start and end: $match_end";
	    }
	}
    }
}
close $fas;

