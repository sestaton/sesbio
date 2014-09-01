#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;
use Bio::TreeIO;
use Getopt::Long;

my $infile;
my $outfile;
my $iformat; ## Optional (Default: nexus)
my $oformat; ## Optional (Default: newick)

my $usage = "USAGE: perl $0 -i infile -o outfile -if format -of format\n\n".
    "Options 'if' and 'of' are optional (Defaults: nexus for 'if' and newick for 'of')."; 

GetOptions(
           'i|infile=s'     => \$infile,
           'o|outfile=s'    => \$outfile,
           'if|informat:s'  => \$iformat,
           'of|outformat:s' => \$oformat,
           );

say $usage and exit(0) if !$infile or !$outfile;

$iformat //= 'nexus';
$oformat //= 'newick';

my $treein  = Bio::TreeIO->new(-file => $infile,      -format => $iformat);
my $treeout = Bio::TreeIO->new(-file => ">$outfile" , -format => $oformat);

while ( my $tree = $treein->next_tree() ) {
    $treeout->write_tree($tree);
}
