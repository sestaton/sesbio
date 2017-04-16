#!/usr/bin/env perl

use strict;
use warnings;

my $base = "http://phylotastic-wg.cs.nmsu.edu/script/phylotastic.cgi";

my ($tree, $taxa) = @ARGV;

$taxa =~ s/[ _]+/+/g;
$taxa =~ s/,/%2C/g;

system("curl \"$base?species=$taxa&tree=$tree&format=newick\" > out.tre");
