#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;
use Bio::DB::HTS::Kseq;

my $usage  = "perl $0 infile > outfile\n";
my $infile = shift or die $usage;

my $kseq = Bio::DB::HTS::Kseq->new($infile);
my $iter = $kseq->iterator;

my %families;
while (my $seqio = $iter->next_seq) {
    my $id = $seqio->name;
    my $famname = ($id =~ /(^\w{3}_\w+?)_/) ? $1 : $id;
    #say $famname;
    push @{$families{ $famname }}, $id;
}

for my $fam (reverse sort { @{$families{$a}} <=> @{$families{$b}} } keys %families) {
    my $ct = @{$families{$fam}};
    say join "\t", $fam, $ct;
}

