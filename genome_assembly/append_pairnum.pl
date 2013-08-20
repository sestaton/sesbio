#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;
use autodie qw(open);
use Getopt::Long;

my $usage = "\n$0 -i in_fas -o outfile -n num\n\n";
my $infile;
my $outfile;
my $num;

GetOptions(
	   'i|infile=s'  => \$infile,
	   'o|outfile=s' => \$outfile,
	   'n|num=i'     => \$num,
	   );

die $usage if !$infile or !$outfile or !$num;

open my $in, '<', $infile;
open my $out, '>', $outfile;

my $totalct = 0;
my $pair;

if ($num == 1) {
    $pair = "/1";
}
elsif ($num == 2) {
    $pair = "/2";
}
else {
    die "\nERROR: Could not determine pair information (must be 1 or 2). Exiting.\n";
}

{
    local $/ = '>';

    while (my $line = <$in>) {
	chomp $line;
	my ($seqid, @seqparts) = split /\n/, $line;
	my $seq = join '', @seqparts;
	next unless defined $seqid && defined $seq;
	$totalct++;
	say $out join "\n", ">".$seqid.$pair, $seq;
    }
}

say "\n=====> $totalct total reads were written to $outfile.\n";

