#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;
use autodie qw(open);
use File::Basename;
use Getopt::Long;
use Data::Dump;
use Bio::Tools::GFF;

my %opt;
my %genes;

GetOptions(\%opt, 'infile|i=s', 'fasta|f=s', 'outfile|o=s');

usage() and exit(0) if !$opt{infile} or !$opt{fasta} or !$opt{outfile};

my $gffio = Bio::Tools::GFF->new( -file => $opt{infile}, -gff_version => 3 );

while (my $feature = $gffio->next_feature()) {
    if ($feature->primary_tag eq 'gene') {
	my @string = split /\t/, $feature->gff_string;
	my ($id) = ($string[8] =~ /ID=?\s+?(gene\d+)/);
	my ($start, $end) = ($feature->start, $feature->end);
	$genes{$id} = join "||", $string[0], $start, $end;
    }
}

#dd \%genes and exit;
open my $out, '>', $opt{outfile};

for my $gene (sort keys %genes) {
    my ($src, $start, $end) = split /\|\|/, $genes{$gene};
    my $tmp = $gene.".fasta";
    system("samtools faidx $opt{fasta} $src:$start-$end > $tmp");
    my $id = join "_", $gene, $src, $start, $end;

    if (-s $tmp) {
	my $seqio = Bio::SeqIO->new( -file => $tmp, -format => 'fasta' );
	while (my $seqobj = $seqio->next_seq) {
	    my $seq = $seqobj->seq;
	    if ($seq) {
		$seq =~ s/.{60}\K/\n/g;
		say $out join "\n", ">$id", $seq;
	    }
	}
    }
    unlink $tmp;
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
