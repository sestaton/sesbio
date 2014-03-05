#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;

while (my $line = <>) {
    chomp $line;
    my $log = $line / 10000;
    my $rounded = int($log + $log/abs($log*2));
    say $rounded;
}
