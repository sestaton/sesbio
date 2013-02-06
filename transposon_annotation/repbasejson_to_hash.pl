#!/usr/bin/env perl

use 5.014;
use strict;
use warnings;
use autodie qw(open);
use lib qw(/home/jmblab/statonse/apps/perlmod/Data-Dump-1.21/blib/lib);
use Data::Dump qw(dd);
use JSON;

my $usage = "$0 repbase1801.json\n";
my $infile = shift or die $usage;
open(my $in, '<', $infile);

my $json = JSON->new->utf8->space_after->decode($in);
close($in);
dd $json;
