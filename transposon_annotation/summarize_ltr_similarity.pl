#!/usr/bin/env perl

## simple script to print (tab-delimited) LTR information:
## LTRname   Length    Percent_identity_of_LTRs

use 5.010;
use strict;
use warnings;
use autodie;

my $usage = "$0 tephra_ltrs.gff3";
my $gff   = shift or die $usage;

open my $in, '<', $gff;

while (<$in>) {
    chomp;
    next if /^#/; 
    my @f = split /\t/; 
    if ($f[2] eq 'LTR_retrotransposon') { 
	my $len = $f[4] - $f[3] + 1; 
	my ($id, $pid) =  ($f[8] =~ /Parent=(repeat_region\d+)\;ltr_similarity=(\d+\.\d+)/); 
	say join "\t", $id, $len, $pid;
    }
}
close $in;
