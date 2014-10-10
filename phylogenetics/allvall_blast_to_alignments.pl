#!/usr/bin/env perl

use 5.012;
use strict;
use warnings;
use open qw(:encoding(UTF-8) :std);
use autodie qw(open);
use utf8;
use charnames qw(:full :short);
use Getopt::Long;
use Data::Dump qw(dd);

my $infile;  # m8 blast
my $dir;     # dir containing the fasta files
my $usage = "perl $0 -i allvall.bln -d dirname\n";

GetOptions(
           'i|infile=s'  => \$infile,
           'd|dir:s'     => \$dir,
           );

die $usage if !$infile; # or !$dir;

my %list;
my %seen;
my %seenall;
my %matchrange;

my @fas_files = glob "./*.fasta";
if (scalar @fas_files < 1) {
    say "\nERROR: Could not find any fasta files in current directory. Exiting.\n";
    exit(1);
}

open my $in, '<', $infile;
while (<$in>) {
    chomp;
    my @f = split;
    if (@f) {
	next if $f[0] eq $f[1];
	next if $f[3] < 2000;
	my $k = mk_key($f[0], $f[1]);
	next if exists $seen{$k};
	#next if exists $seenall{$f[1]};
	if (exists $list{$f[0]}) {
	    push @{$list{$f[0]}}, mk_key($f[1], $f[6], $f[7], $f[8], $f[9]);
	}
	else {
	    $list{$f[0]} = [ mk_key($f[1], $f[6], $f[7], $f[8], $f[9]) ];
	}
	$seen{$k} = 1;
	#$seenall{$f[0]} = 1;
	#$seenall{$f[1]} = 1;
    }
}
close $in;
#dd \%list;

for my $k (keys %list) {
    if (scalar @{$list{$k}} == 16) {
	my @syn_block;
	#my @syn_block = @{$list{$k}}; push @syn_block, $k;
	#my @syn_block_sort = sort { $a cmp $b } @syn_block;
	for my $v (@{$list{$k}}) {
	    my ($t, $qs, $qe, $ts, $te) = mk_vec($v);
	    say join "\t", $k, $t, $qs, $qe, $ts, $te;
	}
	say "#######";
    }
}

#
# methods
#
sub mk_key { join "\N{INVISIBLE SEPARATOR}", map { $_ // " " } @_ }

sub mk_vec { split "\N{INVISIBLE SEPARATOR}", shift }
