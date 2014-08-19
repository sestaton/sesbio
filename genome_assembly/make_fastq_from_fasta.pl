#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;
use autodie qw(open);
use Getopt::Long;

my $usage = "\n$0 -i in_fas -o outfile\n\n";
my $infile;
my $outfile;

GetOptions(
	   'i|infile=s'  => \$infile,
	   'o|outfile=s' => \$outfile,
	   );

die $usage if !$infile or !$outfile;

open my $in, '<', $infile;
open my $out, '>', $outfile;

{
    local $/ = '>';

    while (my $line = <$in>) {
	chomp $line;
	my ($seqid, @seqparts) = split /\n/, $line;
	my $seq = join '', @seqparts;
	next unless defined $seqid && defined $seq;
	say $out join "\n", "@".$seqid, $seq, "+", "I" x length($seq);
    }
}

