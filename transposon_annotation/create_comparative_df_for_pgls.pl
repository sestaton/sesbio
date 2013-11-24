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
    #say "$spec => $bp";
    $cvalh->{$spec} = $bp;
}

for my $file (@files) {
    my ($species) = split(/\_/, $file, 2);
    open my $in, '<', $file;
    while (<$in>) {
	chomp;
	next if /^ReadNum/;
	my @f = split "\t";
	$fams{$f[2]} = 1;
	my $key = mk_key($f[2], $species);
	$sph{$key} = 1;
	#$sph{$f[1]}{$species} = 1;
	#$df{$species}{$f[1]}{$f[2]} = $f[5];
	$df{$species}{$f[2]} = $f[5];
    }
    close $in;
}

#dd \%df; exit;
#my %sfdf;
my $fam_ct = scalar keys %fams;
print "Species\t", join "\t", (sort keys %fams), "\n";
keys %fams;

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

##TODO: Fix bug with creating format that can't be read in R without reformatting in excel
for my $species (sort keys %df) {
    print $species_map{$species};
    for my $f (sort keys %fams) {
	my $key = mk_key($f, $species);
	if (exists $sph{$key} && exists $df{$species}{$f} && exists $cvalh->{$species}) {
	    my $bpsize = $df{$species}{$f} * $cvalh->{$species};
	    #my $log_rounded = int(log($bpsize) + log($bpsize)/abs(log($bpsize)*2));
	    #print "\t$log_rounded";
	    my $rounded = int($bpsize + $bpsize/abs($bpsize*2));
	    print "\t$rounded";
	}
	else {
	    print "\t0";
	}
    }
    print "\n";
}
#dd \%sfdf;
#print $fam_ct," expected families\n";

## subs
sub mk_key { join "\N{INVISIBLE SEPARATOR}", map { $_ // " " } @_ }

sub mk_vec { split "\N{INVISIBLE SEPARATOR}", shift }
