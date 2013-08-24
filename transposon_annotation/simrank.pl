#!/usr/bin/env perl

## TODO: Store matches with Storable for later
use 5.010;
use strict;
use warnings;
use String::Simrank;
use Data::Dump qw(dd dump);
use Getopt::Long;
use autodie qw(open);
use Cwd;

my $usage = "$0 -d infile(db) -q query -o outfile\n";
my $infile; # = shift or die $usage;
my $query; # = shift or die $usage;
my $outfile; # = shift or die $usage;
my $cwd = getcwd();

GetOptions(
           'q|query=s'   => \$query,
           'd|db=s'      => \$infile,
           'o|outfile=s' => \$outfile,
          );

die $usage if !$infile || !$outfile || !$query;

open my $out, '>', $outfile;

my $sr = new String::Simrank( { data => $infile} );
$sr->formatdb({wordlen => 8, valid_chars => 'ACGT', silent => 1});
my $matches = $sr->match_oligos({ query => $query,
				  outlen => 10,
				  minpct => 50,
				  reverse => 1,});

#print $out dump($matches); ## for debug

for my $k (sort { scalar(@{$matches->{$a}}) <=> scalar(@{$matches->{$b}}) } keys %$matches) {
    my $num_matches = scalar @{$matches->{$k}};
    if ($num_matches > 1) { 
	say {$out} join "\t", "$k:", $num_matches;
	for my $hit (@{$matches->{$k}}) { ## Matches are already sorted by percent identity, so not sorting necessary.
	    say {$out} join "\t", "HitID:", $hit->[0], "Perc:", $hit->[1];
	}
    }
}
close $out;
