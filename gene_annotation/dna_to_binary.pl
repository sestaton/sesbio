#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;

my $str = q{ATCGTAGCTACGATCGT};
my $bin = unpack('B*', $str);

my $char = length($bin);
my @arr =  pack("B$char", $bin);

say join "\n", "Orig:         $str", "Binary:       $bin", "Orig (again): @arr";

#Orig:         ATCGTAGCTACGATCGT
#Binary:       0100000101010100010000110100011101010100010000010100011101000011010101000100000101000011010001110100000101010100010000110100011101010100
#Orig (again): ATCGTAGCTACGATCGT
