#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;
use autodie qw(open);
use Getopt::Long;
use Data::Dump;

my %opt;
my %genes;

GetOptions(\%opt, 'infile|i=s', 'fasta|f=s');

my $usage = "$0 -i file.gff -f seq.fasta\n";
die $usage if !$opt{infile} or !$opt{fasta};

open my $in, '<', $opt{infile};
while (<$in>) {
    chomp;
    next if /^#/;
    my @f = split;
    my $id;
    if ($f[2] eq 'gene') {
	($id = $f[8]) =~ /ID=(\w+);/; 
	$id =~ s/ID=//;
	$genes{$id} = join "-", $f[3], $f[4];
    }
}
close $in;

#dd \%genes;
my $seq = seq_to_str($opt{fasta});

for my $gene (keys %genes) {
    my ($start, $end) = split /\-/, $genes{$gene};
    my $length = $end - $start;

    my $outfile = $gene.".fasta";
    open my $out, '>', $outfile;
    my $gene_seq = substr $seq, $start, $length;
    $gene_seq =~ s/.{60}\K/\n/g;
    say join ", ", $gene, $start, $end, length($gene_seq);
    say $out join "\n", ">".$gene."_".$start."-".$end, $gene_seq;
}

sub seq_to_str {
    my ($fasta) = @_;
    open my $fas, '<', $fasta;
    my $seq;

    while (my $line = <$fas>) {
	chomp $line;
	unless ($line =~ /^>/) {
	    $seq .= $line;
	}
    }
    close $fas;

    return $seq;
}
