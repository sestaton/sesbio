#!/usr/bin/env perl

use 5.020;
use warnings;
use autodie;
use File::Basename;
use Data::Printer;

my $usage = "$0 clusterfile\n";
my $file  = shift or die $usage;

my (%cls, $clusnum);
open my $in, '<', $file;
while (my $line = <$in>) {
    chomp $line;
    next if $line =~ /^#/;
    if ($line =~ /^(\d+):/) {
	$clusnum = $1;
    }
    elsif ($line =~ /^\s+(\S+)/) {
	push @{$cls{$clusnum}}, $1;
    }
}

p %cls;
