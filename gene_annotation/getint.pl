## this is a dead simple demonstration of how to get
## overlaps for an interval range

use 5.010;
use strict;
use warnings;
use Set::IntervalTree;
use Data::Dump qw(dd);

my $tree = Set::IntervalTree->new;
while (<DATA>) {
    my @ints = split;
    $tree->insert($ints[1],$ints[2],$ints[3]);
}

my $res = $tree->fetch(10, 30);
say "Found: ",scalar(@$res)," overlaps.";
dd $res;

__DATA__
chr1    gene1    25    30
chr1    gene2    15    20
chr1    gene3    80    90
chr1    gene4    34    37
chr1    gene5    25    28
chr1    gene6    10    13
