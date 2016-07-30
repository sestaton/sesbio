#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;
use Bio::DB::HTS::Faidx;

my $usage = "$0 seq";
my $seq = shift or die $usage;

my $id = 'Zea_mays_chr4_trna181-AlaAGC';
my $index = Bio::DB::HTS::Faidx->new($seq);

my ($sseq, $len) = $index->get_sequence($id);
say join "\n", ">$id", $sseq;
