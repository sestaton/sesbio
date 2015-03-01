#!/usr/bin/env perl

use 5.020;
use strict;
use warnings;
use autodie;
use File::Find;
use File::Basename;
use Data::Dump;
use List::MoreUtils qw(natatime);
use experimental 'signatures';

my $usage = "$0 dir_of_reports\n";
my $dir   = shift or die $usage;

my @files;
find( sub { push @files, $File::Find::name if -f and /\.txt_summary.html$/ }, $dir );

my @calls = grep { /calls/ } @files;
my @filt  = grep { /filt/  } @files;
#dd \@calls; # and exit;
#dd \@filt and exit;

my %allstats;
for my $file (@calls) {
    my ($acc) = ($file =~ /^$dir\/(\w+)_/);
    my $stats = parse_vep($file);
    $allstats{$acc} = $stats;
}

#dd \%allstats;
for my $samp (keys %allstats) {
    for my $tab (keys %{$allstats{$samp}}) {
	my $outfile = $dir."_".$tab."_calls_vep_summary.tsv";
	open my $out, '>>', $outfile;
	#say $out join "\t", "Accession", "Type", "Count";
	for my $ty (@{$allstats{$samp}{$tab}}) {
	    for my $k (keys %$ty) {
		say $out join "\t", $samp, $k, $ty->{$k};
	    }
	}
	close $out;
    }
}

sub parse_vep ($file) {
    my %stats;

    open my $in, '<', $file;
    while (my $l = <$in>) {
	chomp $l;
	next unless $l =~ /\S/;
	if ($l =~ /var (\w+) = drawTable/) {
	    my $type = $1;
	    next if $type eq 'chr_table';
	    my @res = ($l =~ /\[\'(\w+)\'\,(\d+)\]/g);
	    next unless @res;
	    my $it = natatime 2, @res;
	    while (my @vals = $it->()) {
		my ($k, $v) = @vals;
		push @{$stats{$type}}, { $k => $v };
	    }
	}
    }
    close $in;

    return \%stats;
}
