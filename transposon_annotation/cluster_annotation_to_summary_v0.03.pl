#!/usr/bin/env perl

use v5.12;
use strict;
use warnings;
use autodie qw(open);
use Data::Dump qw(dd dump);
use List::Util qw(sum);

my $usage = "$0 annot.tsv\n";
my $infile = shift or die $usage;

my %annot;
my %sfam_hitct;
my %fam_readct;
my $total_ct = 0;

open(my $in, '<', $infile);

while (<$in>) {
    chomp;
    next if /^Cluster/;
    my @fields = split;
    $total_ct += $fields[1];
    if (scalar @fields == 4) {
	$sfam_hitct{$fields[2]}++;
	$fam_readct{$fields[3]} += $fields[1];
	if (exists $annot{$fields[2]}{$fields[3]}) {
	    #$sfam_hitct{$fields[2]}++;
	    #$fam_readct{$fields[3]} += $fields[1];
	    push @{$annot{$fields[2]}{$fields[3]}}, $fields[1];
	}
	else {
	    #$sfam_hitct{$fields[2]} = 1;
	    #$fam_readct{$fields[3]} = $fields[1];
	    $annot{$fields[2]}{$fields[3]} = [$fields[1]];
	}
    }
    else {
	my $fam = $fields[5];
	$fam =~ s/\-\d.*//;
	$sfam_hitct{$fields[4]}++;
	$fam_readct{$fam} += $fields[1];
	if (exists $annot{$fields[4]}{$fam}) {
	    #$sfam_hitct{$fields[4]}++;
            #$fam_readct{$fam} += $fields[1];
	    push @{$annot{$fields[4]}{$fam}}, $fields[1];
	}
	else {
	    #$sfam_hitct{$fields[4]} = 1;
            #$fam_readct{$fam} = $fields[1];
	    $annot{$fields[4]}{$fam} = [$fields[1]];
	}
    }
}
close($in);

#dd \%annot;
#dd \%sfam_hitct;
#dd \%fam_readct;
for my $sf (reverse sort { $sfam_hitct{$a} <=> $sfam_hitct{$b} } keys %sfam_hitct) {
    for my $f (reverse sort { $fam_readct{$a} <=> $fam_readct{$b} } keys %fam_readct) {
	if (exists $annot{$sf}{$f}) {
	    my $read_ct = sum @{$annot{$sf}{$f}};
	    my $perc_cov = sprintf("%.12f",$read_ct/$total_ct);
	    say join "\t", $sf, $f, $read_ct."/".$total_ct, $perc_cov;
	}
    }
}

