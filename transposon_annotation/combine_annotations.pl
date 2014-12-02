#!/usr/bin/env perl

## This script combines annotations from Transposome, and this
## may be useful for looking at variation from different random
## read samples at the same level of coverage.

use 5.010;
use strict;
use warnings;
use autodie qw(open);
use Statistics::Descriptive;
use Getopt::Long;

my $usage = "$0 -i file file file\n";
my @files;
my %res;

GetOptions( 'i|infiles=s{1,}' => \@files );

if (@files < 1) {
    say "\nERROR: Command line not parsed correctly. Incorrect number of parameters passed. Exiting.\n";
    say $usage;
    exit(1);
}

for my $file (@files) {
    open my $in, '<', $file;
    while (<$in>) {
	chomp;
	next if /^ReadNum/;
	#ReadNum Superfamily Family ReadCt/ReadsWithHit HitPerc GenomePerc
	my @f = split;
	if (exists $res{$f[1]}{$f[2]}) {
	    push @{$res{$f[1]}{$f[2]}}, $f[5];
	}
	else {
	    $res{$f[1]}{$f[2]} = [ $f[5] ];
	}
    }
    close $in;
}

say join "\t", "Superfamily", "Family", "Mean", "SD";

for my $sfam (sort keys %res) {
    for my $fam (sort keys %{$res{$sfam}}) {
	if (scalar @{$res{$sfam}{$fam}} > 1) {
	    ## add data
	    my $stat = Statistics::Descriptive::Full->new;
	    $stat->add_data(@{$res{$sfam}{$fam}});
	    ## calc mean, sd
	    my $mean = $stat->mean;
	    my $sd   = $stat->standard_deviation;
	    my $var  = $stat->variance;
	    ## print results
	    say join "\t", $sfam, $fam, $mean, $sd;
	    undef $stat;
	}
	else {
	    # not all superfamilies are represented more than once, so the std dev is 0
	    say join "\t", $sfam, $fam, $res{$sfam}{$fam}->[0], 0;
	}
    }
}
	
