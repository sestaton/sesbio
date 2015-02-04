#!/usr/bin/env perl

##NB: This takes the first N bases of a sequence

use 5.010;
use strict;
use warnings;
use autodie qw(open);
use Getopt::Long;

my $usage = "\nUSAGE: perl $0 -i in_fas -o outfile -l length\n\n";
my $infile;
my $outfile;
my $length;

GetOptions(
	   'i|infile=s'  => \$infile,
	   'o|outfile=s' => \$outfile,
	   'l|length=i'  => \$length,
	   );

die $usage if !$infile or !$outfile or !$length;

open my $in, '<', $infile;
open my $out, '>', $outfile;

my ($totalct, $subseqct) = (0, 0);

{
    local $/ = '>';

    while (my $line = <$in>) {
	chomp $line;
	$totalct++;
	my ($seqid, @seqparts) = split /\n/, $line;
	my $seq = join '', @seqparts;
	next unless defined $seqid && defined $seq;
	next unless length($seq) >= $length;
	my $subseq = substr $seq, 0, $length;
	$subseqct++;
	say $out join "\n", ">".$seqid, $subseq;
    }
}
close $in;
close $out;

say "\n=====> $subseqct of $totalct total reads were written to $outfile.\n";

