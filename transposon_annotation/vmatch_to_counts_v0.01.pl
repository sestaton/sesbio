#!/usr/bin/env perl

use 5.014;
use strict;
use warnings;
use autodie qw(open);
use lib qw(/home/jmblab/statonse/apps/perlmod/Data-Dump-1.21/blib/lib);
use Data::Dump qw(dd);
use JSON;
use List::Util qw(sum);

my $usage = "$0 idlist vmatches\n";
my $infile = shift or die $usage;
my $vmatches = shift or die $usage;
open(my $in, '<', $infile);
open(my $vm, '<', $vmatches);

my %family_map;

while (my $line = <$in>) {
    chomp $line;
    my ($f, $sf, $source)  = split /\t/, $line;
    next unless defined $sf && defined $f; ## why?
    if ($sf =~ /(\s+)/) {
	$sf =~ s/$1/\_/;
    }
    $f =~ s/\s/\_/;
    if (exists $family_map{$sf}) {
	push @{$family_map{$sf}}, {$f => 0};
    }
    else {
	$family_map{$sf} = [];
    }
}
close($in);

while (my $l = <$vm>) {
    chomp $l;
    $l =~ s/^\s+//;
    my ($ct, $fam_match)  = split /\s/, $l;
    #say join("\t", $ct, $fam_match);
    if ($ct > 1) {
	for my $superfamily (keys %family_map) {
	    #say $superfamily;
	    while (my ($family_index, $family) = each @{$family_map{$superfamily}}) {
		#say $family;
		for my $fam_name (keys %$family) {
		    #say $fam_name;
		    if (exists $family->{$fam_match}) {
			#say "$fam_name and $fam_match match";
			#push @{$family_map{$superfamily}[$family_index]{$fam_name}}, $ct;
			$family_map{$superfamily}[$family_index]{$fam_name} = $ct;
		    }
		}
	    }
	}
    }
}
close($vm);

my %supfamily_match_ct;

for my $supfam (keys %family_map) {
    while (my ($fam_index, $fam) = each @{$family_map{$supfam}}) {
	for my $fam_n (keys %$fam) {
	    if ($fam->{$fam_n} > 0) {
		#say $supfam,' ',$fam_n,' ',$fam->{$fam_n};
		if (exists $supfamily_match_ct{$supfam}) {
		    push @{$supfamily_match_ct{$supfam}}, $fam->{$fam_n};
		#    say $supfam,' is in da hash!';
		}
		else {
		    $supfamily_match_ct{$supfam} = [];
		#    #say $supfam, ' is not in da hash';
		}
	    }
	}
    }
    my $sum = sum(@{$supfamily_match_ct{$supfam}});
    if (defined $sum) {
	say join("\t",$supfam,$sum);
    }
}


#for my $key (keys %supfamily_match_ct) {
#    my $sum = sum(@{$supfamily_match_ct{$key}});    

				

#dd %supfamily_match_ct;
