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
    if ($line =~ /^(\d+):/) { # cluster number
	$clusnum = $1;
    }
    elsif ($line =~ /^\s+(\S+)/) { #reads in cluster, one per line
	push @{$cls{$clusnum}}, $1;
    }
}
close $in;

p %cls; # a nice hash with an array of read ids for values
