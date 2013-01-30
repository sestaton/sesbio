#!/usr/bin/env perl

use strict;
use warnings;
use String::Simrank;
use Data::Dump qw(dd dump);
use Getopt::Long;
use feature 'say';
use autodie qw(open);
use Cwd;

my $usage = "$0 -d infile(db) -q query -o outfile\n";
my $infile; # = shift or die $usage;
my $query; # = shift or die $usage;
my $outfile; # = shift or die $usage;
my $cwd = getcwd();
#say "CWD:  $cwd";

GetOptions(
           'q|query=s' => \$query,
           'd|db=s'    => \$infile,
           'o|outfile=s' => \$outfile,
          );

if (!$infile || !$outfile || !$query) { die $usage; }

#my $out = $cwd."/".$outfile;
open(my $out, '>', $outfile);

my $sr = new String::Simrank( { data => $infile} );
$sr->formatdb({wordlen => 8, valid_chars => 'ACGT', silent => 1});
my $matches = $sr->match_oligos({ query => $query,
				  outlen => 10,
				  minpct => 50,
				  reverse => 1,});
				  #outfile => $out,});

print $out dump($matches);

#for my $k (keys %$matches) {
#for my $k (reverse sort { $matches->{$a}[0] <=> $matches->{$b}[0] } keys %$matches) {
#    my $num_matches = scalar @{$matches->{$k}};
#    if ($num_matches > 1) { # number of matches
#	print {$out} "$k:\t$num_matches\t";
#	for my $hit (@{$matches->{$k}}) {
#	    print {$out} "HitID:\t".$hit->[0]."\tPerc:\t".$hit->[1],"\n";
#	}
#    }
#}
close($out);
