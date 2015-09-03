#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;
use File::Basename;
use autodie qw(open);
use Getopt::Long;
use List::Util qw(min max);
use Data::Dump;

my %opt;
my %genes;

my $script = basename($0, ());
GetOptions(\%opt, 'infile|i=s', 'fasta|f=s', 'debug');

my $usage = "$script -i file.gff -f seq.fasta\n";
die $usage if !$opt{infile} or !$opt{fasta};

my ($id, $eid);
my $exons = 0;
open my $in, '<', $opt{infile};
while (<$in>) {
    chomp;
    next if /^#/;
    my @f = split /\t/;

    if ($f[2] eq 'gene') {
	$exons = 0;
	($id = $f[8]) =~ /ID=(\w+);/; 
	$id =~ s/ID=//;
    }

    if ($f[2] eq 'exon') {
	$exons++;
	($eid = $f[8]) =~ /Parent=(\w+)/;
        $eid =~ s/ID=.*;//;
	$eid =~ s/Parent=//;
	if ($id eq $eid) {
	    my $exonid = "exon".$exons;
	    $genes{$id}{$exonid} = join "-", $f[3], $f[4];
	}
    }
}
close $in;

my @range;
dd \%genes if $opt{debug};
my $seq = seq_to_str($opt{fasta});

for my $gene (keys %genes) {
    my $outfile = $gene."_CDS.fasta";
    open my $out, '>>', $outfile unless $opt{debug};
    my $cds;

    my @exons = map  { $_->[0] }
     	        sort { $a->[1] <=> $b->[1] }
                map  { [ $_, /(\d+)/ ] }
                keys %{$genes{$gene}};

    for my $exon (@exons) {
	my ($start, $end) = split /\-/, $genes{$gene}{$exon};
	say join ", ", $gene, $exon, $start, $end if $opt{debug};
	push @range, $start, $end;
	my $length = $end - $start;
	my $cds_seq = substr $seq, $start, $length;
	$cds .= $cds_seq;
    }
    my $min = min(@range);
    my $max = max(@range);
    $cds =~ s/.{60}\K/\n/g;
    say join ", ", $gene, $min, $max, length($cds) if $opt{debug};
    say $out join "\n", ">".$gene."_".$min."-".$max, $cds unless $opt{debug};

    close $out unless $opt{debug};
    undef $cds;
    undef @range;
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
