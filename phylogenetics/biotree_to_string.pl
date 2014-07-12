#!/usr/bin/env perl

## This is my answer to this question on SO: http://stackoverflow.com/a/24655606/1543853

use strict;
use warnings;
use Bio::TreeIO;

my $treeio  = Bio::TreeIO->new(-format => 'newick', -fh => \*DATA);
my $treeout = Bio::TreeIO->new(-format => 'newick', -fh => \*STDOUT);

while (my $tree = $treeio->next_tree) {
    $treeout->write_tree($tree);
}

__DATA__
(A:9.70,(B:8.234,(C:7.932,(D:6.321,((E:2.342,F:2.321):4.231,((((G:4.561,H:3.721):3.9623,
I:3.645):2.341,J:4.893):4.671)):0.234):0.567):0.673):0.456);
