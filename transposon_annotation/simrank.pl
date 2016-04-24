#!/usr/bin/env perl

## The module required by this script: https://metacpan.org/pod/release/SHURIKO/String-Simrank-0.079/lib/String/Simrank.pm
## Reference: http://bmcecol.biomedcentral.com/articles/10.1186/1472-6785-11-11

use 5.010;
use strict;
use warnings;
use autodie;
use String::Simrank;
use Getopt::Long;
use Cwd;

my $usage = "$0 -d infile(db) -q query -o outfile\n";
my $infile;
my $query;
my $outfile;
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

for my $k (sort { scalar(@{$matches->{$a}}) <=> scalar(@{$matches->{$b}}) } keys %$matches) {
    my $num_matches = scalar @{$matches->{$k}};
    if ($num_matches > 1) { 
	say {$out} join "\t", "$k:", $num_matches;
	## Matches are already sorted by percent identity, so no sorting necessary.
	for my $hit (@{$matches->{$k}}) {
	    say {$out} join "\t", "HitID:", $hit->[0], "Perc:", $hit->[1];
	}
    }
}
close $out;
