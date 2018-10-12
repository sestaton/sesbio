#!/usr/bin/env perl

#TODO: Handle compressed input 

use 5.010;
use strict;
use warnings;
use Scalar::Util qw(openhandle);
use File::Basename;
use File::Spec;
use Sort::Naturally;

my $usage = "$0 genome.gff\n";
my $gfffile = shift or die $usage;

split_gff($gfffile);

exit;
#
# methods
#
sub split_gff {
    my ($gff) = @_;

    my ($name, $path, $suffix) = fileparse($gff, qr/\.[^.]*/);
    open my $in, '<', $gff or die "\nERROR: Could not open file: $gff\n";
    #open my $in, '-|', 'zcat', $gff or die $!;

    my $header;
    my %regions;
    while (my $line = <$in>) {
	chomp $line;
	next if $line =~ /^###$/;
	if ($line =~ /^##/) {
	    next if $line =~ /gff\-version/;
	    my ($reg, $id, $start, $end) = split /\s+/, $line;
	    $regions{$id} = $end;
	}
	else {
	    last;
	}
    }
    close $in;

    open my $gffio, '<', $gff or die "\nERROR: Could not open file: $gff\n";
    #open my $gffio, '-|', 'zcat', $gff or die $!;

    my ($outfile, $out, %seen);
    while (my $line = <$gffio>) {
        chomp $line;
        next if $line =~ /^#/;
	my @f = split /\t/, $line;
	unless (exists $seen{$f[0]}) {
	    $outfile = File::Spec->catfile($path, $name.'_'.$f[0].'.gff3');
	    close $out if openhandle($out);
	    open $out, '>', $outfile or die $!;
	    say $out '##gff-version 3';
	    say $out join q{ }, '##sequence-region', '1', $f[0], $regions{$f[0]};
	    $seen{$f[0]} = 1;
	}
	say $out join "\t", @f;
    }
    close $gffio;

    return;
}
