#!/usr/bin/env perl

use strict;
use warnings;
use Bio::TreeIO;

my $usage = "USAGE: perl $0 infile outfile\n"; 
my $infile = $ARGV[0] or die $usage;
my $outfile = $ARGV[1] or die $usage;

my $treein  = Bio::TreeIO->new(-file => $infile, '-format' => 'nexus');
my $treeout = Bio::TreeIO->new(-file => ">$outfile" , '-format' => 'newick');

while ( my $tree = $treein->next_tree() ) {
    $treeout->write_tree($tree);
}
