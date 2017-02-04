#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;
use autodie qw(open);
use File::Basename;
use Bio::DB::HTS::Faidx;
use Sort::Naturally;
use Getopt::Long;

my %opt;
my %genes;

GetOptions(\%opt, 'infile|i=s', 'fasta|f=s', 'outfile|o=s');

usage() and exit(0) if !$opt{infile} or !$opt{fasta} or !$opt{outfile};

open my $in, '<', $opt{infile};

while (my $line = <$in>) {
    chomp $line;
    next if $line =~ /^#/;
    last if $line =~ />/;
    my @string = split /\t/, $line; 
    next if $string[0] =~ /HanXRQCP|HanXRQMT/;
    die $line unless defined $string[2];
    if ($string[2] eq 'gene') {
	my ($id) = ($string[8] =~ /ID=(\w+\d+);/);
	die $line unless defined $id;
	my ($start, $end) = ($string[3], $string[4]);
	$genes{$id} = join "||", $string[0], $start, $end;
    }
}
close $in;

my $index = Bio::DB::HTS::Faidx->new($opt{fasta});

open my $out, '>', $opt{outfile};

for my $gene (nsort keys %genes) {
    my ($src, $start, $end) = split /\|\|/, $genes{$gene};
    my $id = join q{ }, $gene, join "_", $src, "$start-$end";

    my $location = "$src:$start-$end";
    my ($seq, $length) = $index->get_sequence($location);
    
    $seq =~ s/.{60}\K/\n/g;
    say $out join "\n", ">".$id, $seq;
}
close $out;

sub usage {
    my $script = basename($0);
  print STDERR <<END

USAGE: $script -i file.gff -f seqs.fas

Required:
 -i|infile    :    GFF file to extract gene coordinates from
 -f|fasta     :    FASTA file to pull the gene regions from.
 -o|outfile   :    The file to place gene sequences.
    
Options:
 -h|help      :    Print usage statement (not implemented).
 -m|man       :    Print full documentation (not implemented).

END
}
