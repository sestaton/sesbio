#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "blastn -outfmt 6 ... | $0 - > blast_msp_recon.txt\n-OR-\n".
    "$0 blasttable.bln > blast_msp_recon.txt\n";
die $usage if !@ARGV;

while (<>) {
    chomp;
    my @f = split;
    die "This blast report is not formatted correctly. Exiting.\n"
	unless @f == 12;
    next if $f[0] eq $f[1];
    printf("%06d %03d %05d %05d %s %05d %05d %s \n", 
	   $f[11], $f[2], $f[6], $f[7], $f[0], $f[8], $f[9], $f[1]);
}
