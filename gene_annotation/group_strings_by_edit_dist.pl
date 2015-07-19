#!/usr/bin/env perl

##TODO: Take patterns from command line

use 5.010;
use strict;
use warnings;
use autodie qw(open);
use String::Approx 'amatch';

my $usage = "$0 infile\n";
my $infile = shift or die $usage;

open my $in, '<', $infile;

my @scores;
my %matches;
my $ct = 0;

while (<$in>) {
    chomp;
    $ct++;
    #$matches{$_}++;
    push @scores, $_;
}

## look for duplicates first 
#for my $k (reverse sort { $matches{$a} <=> $matches{$b} } keys %matches) {
#    if ($k > 1) {
#        say join "\t", $k, $matches{$k};
#    }
#}

for my $i (1..5) {
    say "===== Groupings with an edit distance of $i =====";
    for my $pat (qw(BACBEDEAEACECEEC EBCDEECCABAACEAE EEDAEECCAAAECEDE ACDDEDCCDBEBD CEDAEEDBBABA)) {
	my @match = amatch($pat, [$i], @scores);
	say $pat, "\n";
	say join "\n", @match;
	say "=" x 20;
    }
}
