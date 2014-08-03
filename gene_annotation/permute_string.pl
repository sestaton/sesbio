#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;

my $usage  = "USAGE: perl $0 <string>\n";
my $string = shift or die $usage;

# Fischer-Krause ordered permutation generator
sub permute (&@) {
    my $code = shift;
    my @idx = 0..$#_;
    while ( $code->(@_[@idx]) ) {
	my $p = $#idx;
	--$p while $idx[$p-1] > $idx[$p];
	my $q = $p or return;
	push @idx, reverse splice @idx, $p;
	++$q while $idx[$p-1] > $idx[$q];
	@idx[$p-1,$q]= @idx[$q,$p-1];
    }
}

permute { say @_ } split //, $string;
