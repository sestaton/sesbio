#!/usr/bin/env perl

use strict;
use warnings;

my $string = shift;

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

permute { print "@_\n" } split //, $string;
