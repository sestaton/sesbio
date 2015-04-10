#!/usr/bin/env perl

use 5.020;
use warnings;
use autodie;
use experimental 'signatures';

my $usage = "$0 repeatdb.fas\n";
my $infile = shift or die $usage;
open my $in, '<', $infile;

while (my $line = <$in>) {
    chomp $line;
    if ($line =~ /^>/) {
	$line =~ s/>//;
	my ($f, $sf, $source)  = split /\t/, $line;
	next unless defined $sf && defined $f;
	if ($sf =~ /(\s+)/) {
	    $sf =~ s/$1/\_/;
	}
	$f =~ s/\s/\_/;
	my $family = map_family_name($f);

	say join q{ }, $f, $family;
    }
}
close $in;


sub map_family_name ($family) {
    my $family_name;

    if ($family =~ /(^RL[GCX][_-][a-zA-Z]*\d*?[_-]?[a-zA-Z-]+?\d*?)/) {
        $family_name = $1;
        #$family_name .= $2 if $2;
    }
    elsif ($family =~ /(^D[HT][ACHMT][_-][a-zA-Z]+\d*?)/) {
        $family_name = $1;
        #$family_name .= $2 if $2;
    }
    elsif ($family =~ /(^PPP[_-][a-zA-Z]+\d*?)/) {
        $family_name = $1;
        #$family_name .= $2 if $2;
    }
    elsif ($family =~ /(^R[IS][LT][_-][a-zA-Z]+\d*?)/) {
        $family_name = $1;
        #$family_name .= $2 if $2;
    }
    else {
        $family_name = $family;
    }

    $family_name =~ s/_I// if $family_name =~ /_I_|_I$/;
    $family_name =~ s/_LTR// if $family_name =~ /_LTR_|_LTR$/;

    return $family_name;
}
