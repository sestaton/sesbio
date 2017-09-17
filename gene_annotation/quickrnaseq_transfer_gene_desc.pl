#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;
use autodie;
use List::MoreUtils qw(uniq);
use Data::Dump::Color;

my $usage = "$0 table.tsv.gz tab.html";
my $table = shift or die $usage;
my $html  = shift or die $usage;

open my $tin, '-|', 'zcat', $table;

my %descs;
while (my $line = <$tin>) {
    chomp $line;
    my @f = split /\t/, $line;
    my ($id, $desc) = @f[0,12];
    next unless defined $desc;
    $id =~ s/-RA//;
    push @{$descs{$id}}, $desc;
}
close $tin;

#dd \%descs and exit;

open my $hin, '<', $html;
while (my $line = <$hin>) {
    chomp $line;
    if ($line !~ /^Ha\d+/) {
	say $line;
    }
    if ($line =~ /^Ha\d+/) {
	my @f = split /\t/, $line;
	my @fields = grep { length } @f;
	if (exists $descs{$fields[0]}) {
	    my @desc = uniq @{$descs{$fields[0]}};
	    say join "\t", @fields[0..2], (join q{,}, @desc), @fields[3..$#fields];
	}
	else {
	    say join "\t", @fields[0..2], '',@fields[3..$#fields];
	}
    }
}
close $hin;
