#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;
use autodie;

my $usage   = "$0 outdir outfile\n";
my $outdir  = shift or die $usage;
my $outfile = shift or die $usage;
my @files   = glob "$outdir/*.tsv";

open my $out, '>', $outfile;
say $out join "\t", "Line", "Ref", "Pos", "Cov";

for my $file (@files) {
    next if $file eq $outfile;
    say STDERR "working on $file...";
    my $ref = (split /\_/, $file)[0];
    open my $in, '<', $file;
    while (<$in>) {
	chomp;
	my @f = split /\t/;
	say $out join "\t", $ref, @f;
    }
    close $in;
}
close $out;
