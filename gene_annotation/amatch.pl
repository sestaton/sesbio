#!/usr/bin/env perl

use strict;
use warnings;
use integer;

my $Sigma = 256;                   # Size of alphabet.
my @po2   = map { 1 << $_ } 0..31; # Cache powers of two.
my $debug = 1;                     # For the terminally curious.

sub amatch {
    my $P = shift; # Pattern.
    my $k = shift; # Amount of degree of proximity.

    #print join "\n", $P, $k if $debug;

    my $m = length $P; # Size of pattern.
    # If no degree of proximity specified assume 10% of the pattern size.
    $k = (10 * $m) / 100 + 1 unless defined $k;

    # Convert pattern into a bit mask.
    my @T = (0) x $Sigma;
    for (my $i = 0; $i < $m; $i++) {
	$T[ord(substr($P, $i))] |= $po2[$i];
    }

    if ($debug) {
	for (my $i = 0; $i < $Sigma; $i++) {
	    printf "T[%c] = %s\n",
		$i, unpack("b*", pack("V", $T[$i])) if $T[$i];
	}
    }

    my (@s, @r); # s: current state, r: previous state.
    # Initialize previous states.
    for ($r[0] = 0, my $i = 1; $i <= $k; $i++) {
	$r[$i] = $r[$i-1];
	$r[$i] |= $po2[$i-1];
    }

    if ($debug) {
	for (my $i = 0; $i <= $k; $i++) {
	    print "r[$i] = ", unpack("b*", pack("V", $r[$i])), "\n";
	}
    }

    my $n = length();    # Text size.
    my $mb = $po2[$m-1]; # If this bit is lit, we have a hit.

    for ($s[0] = 0, my $i = 0; $i < $n; $i++) {
	$s[0] <<= 1;
	$s[0] |= 1;
	my $Tc = $T[ord(substr($_, $i))]; # Current character.
	$s[0] &= $Tc;                     # Exact matching.
	print "$i s[0] = ", unpack("b*", pack("V", $s[0])), "\n"
	    if $debug;

	for (my $j = 1; $j <= $k; $j++) { # Approximate matching.
	    $s[$j] = ($r[$j] << 1) & $Tc;
	    $s[$j] |= ($r[$j-1] | $s[$j-1]) << 1;
	    $s[$j] |= $r[$j-1];
	    $s[$j] |= 1;
	    print "$i s[$j] = ", unpack("b*", pack("V", $s[$j])), "\n"
		if $debug;
	}
	return $i > $m ? $i - $m : 0 if $s[$k] & $mb; # Match.
	@r = @s;
    }
    return -1; # Mismatch.
}

my $P = @ARGV ? shift : "perl";
my $k = shift if @ARGV;

while (<STDIN>) {
    print if amatch($P, $k) >= 0;
}
