#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use autodie qw(open);

my $usage = "$0 -i infile -sfo supfam_abund -fo fam_abund\n";
my $infile;
my $sf_outfile;
my $f_outfile;
my $gcov = 0;
my %fam;
my %superfam;

GetOptions(
           'i|infile=s'           => \$infile,
           'sfo|supfam_outfile=s' => \$sf_outfile,
           'fo|fam_outfile=s'     => \$f_outfile, 
           );

if (!$infile || !$sf_outfile || !$f_outfile) {
    print $usage;
    exit(1);
}

open my $in, '<', $infile;
open my $sfout, '>', $sf_outfile;
open my $fout, '>', $f_outfile;

while (<$in>) {
    chomp;
    next if /^Superfamily/;
    my ($sf, $f, $rc, $perc) = split;
    $superfam{$sf} += $perc;
    $fam{$f} += $perc;
}
close $in;
print $sfout "#Superfamily abundance for: $infile\n";
for my $sfk (reverse sort { $superfam{$a} <=> $superfam{$b} } keys %superfam) {
    print $sfout join "\t", $sfk, $superfam{$sfk}, "\n";
}
close $sfout;
print $fout "#Family abundance for: $infile\n";
for my $fk (reverse sort { $fam{$a} <=> $fam{$b} } keys %fam) {
    print $fout join "\t", $fk, $fam{$fk}, "\n";
}
close $fout;
