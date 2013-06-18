#!/usr/bin/env perl

use v5.14;
use strict;
use warnings;
use autodie qw(open);
use Data::Dump qw(dd);

my $usage = "perl $0 outfile\n";
my $outfile = shift or die $usage;

my @files = glob("*summary_edit.tsv");

die "\nERROR: Could not get files: $!"
    unless scalar @files > 0;

my %df;
my %sfams;

for my $file (@files) {
    my ($species) = split(/\_/, $file, 2);
    open my $in, '<', $file;
    while (<$in>) {
	chomp;
	next if /^ReadNum/;
	my @f = split "\t";
	$sfams{$f[1]} = 1;
	$df{$species}{$f[1]}{$f[2]} = $f[5];
    }
    close $in;
}

#open my $out, '>', $outfile;

#dd \%df;
my %sfdf;

for my $species (sort keys %df) {
    #print "$species\t";
    for my $superfam (sort keys %{$df{$species}}) {
	#print "$superfam\t";
	for my $fam (keys %{$df{$species}{$superfam}}){
	#    print "$fam($df{$species}{$superfam}{$fam})"
	    $sfdf{$species}{$superfam} += $df{$species}{$superfam}{$fam};
	}
    }
    #print "\n";
}

#print "Species";
#for my $k (sort keys %sfams) {
#    print "\t$k";
#}
#print "\n";

for my $species (sort keys %sfdf) {
    print "$species";
    for my $sf (sort keys %{$sfdf{$species}}) {
	print "\t$sf($sfdf{$species}{$sf})"; 
    }
    print "\n";
}
dd \%sfdf;
