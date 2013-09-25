#!/usr/bin/env perl


## TODO: Change names for printing out full Genus name in plot

use 5.014;
use utf8;
use strict;
use warnings;
use warnings FATAL => "utf8";
use autodie qw(open);
use Data::Dump qw(dd);
use charnames qw(:full :short);

my $usage = "perl $0 outfile\n";
my $outfile = shift or die $usage;

my @files = glob("*summary_edit.tsv");

die "\nERROR: Could not get files: $!"
    unless scalar @files > 0;

my %df;
my %fams;
my %sph;

my $cvalh = {
    Ageratina => '1973',
    Ann1238   => '3600',
    CP        => '1702',
    Calyc     => '1887',
    Dasy      => '2533',
    Gerb      => '2494',
    Gnaph     => '1249',
    Saff      => '1364',
    Senecio   => '1540',
    TKS       => '1871',
};

my $mil = 1_000_000;

for my $spec (keys %$cvalh) {
    my $bp = $cvalh->{$spec} * $mil;
    $cvalh->{$spec} = $bp;
}

for my $file (@files) {
    my ($species) = split(/\_/, $file, 2);
    open my $in, '<', $file;
    while (<$in>) {
	chomp;
	next if /^ReadNum/;
	my @f = split "\t";
	if (exists $df{$species}) {
	    push @{$df{$species}}, { $f[2] => $f[5] };
	}
	else {
	    $df{$species} = [ { $f[2] => $f[5] } ];
	}
    }
    close $in;
}

my %species_map = ( 
    Ageratina => 'Ageratina',
    Ann1238   => 'Helianthus',
    CP        => 'Centrapallus',
    Calyc     => 'Nasanthus',
    Dasy      => 'Fulcaldea',
    Gerb      => 'Gerbera',
    Gnaph     => 'Pseudognaphalium',
    Saff      => 'Carthamus',
    Senecio   => 'Senecio',
    TKS       => 'Taraxacum',
    );


my @fams;

#dd \%df;
#exit;

my %stats;
say join "\t", "Species", "Family","PercCov", "BPCov";
for my $sp (sort keys %df) {
    my $sp_fam_ct = scalar @{$df{$sp}};
    for my $fh (@{$df{$sp}}) {
	#for (my ($f, $bp) = each %$fh) {
	for my $f (sort keys %$fh) {
	    my $bpsize = $fh->{$f} * $cvalh->{$sp};
	    #my $min = min(@{$df{$sp}});
	    #my $max = max(@{$df{$sp}});
	    #my $mean = mean(@{$df{$sp}});
	    #my $median = median(@{$df{$sp}});
	    #$stats{$sp} = join "|", $sp_fam_ct, $min, $max, $mean, $median;
	    #$stats{$sp} = join "|", $sp_fam_ct, sprintf("%.8f", $min * 100), sprintf("%.8f", $max * 100), $mean, sprintf("%.8f", $median);
	    say join "\t", $species_map{$sp}, $f, $fh->{$f}, $bpsize;
	}
    }
}

#dd \%stats;

sub min {
    my $min = shift;
    for ( @_ ) { $min = $_ if $_ < $min }
    return $min;
}

sub max {
    my $max = shift;
    for ( @_ ) { $max = $_ if $_ > $max }
    return $max;
}

sub mean { 
    my @array = @_; 
    my $sum; 
    my $count = scalar @array; 
    for (@array) { $sum += $_; } 
    return sprintf("%.8f",$sum / $count); 
}

sub median {
    my @orig_array = @_;
    my @array = sort {$a <=> $b} @orig_array;
    if ($#array % 2 == 0) {
        my $median = $array[($#array / 2)];
        return $median;
    }
    else {
        my $median = $array[int($#array / 2)] + (($array[int($#array / 2) + 1] - $array[int($#array / 2)]) / 2);
        return $median;
    }
}
