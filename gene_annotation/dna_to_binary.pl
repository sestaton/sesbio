#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;

my $str = q{ATCGTAGCTACGATCGT};
my $bin = unpack('B*', $str);

my $char = length($bin);
my @arr =  pack("B$char", $bin);

say join "\n", "Orig:         $str", "Binary:       $bin", "Orig (again): @arr";


