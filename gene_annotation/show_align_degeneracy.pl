#!/usr/bin/env perl

## this was my response to a biostars question: http://www.biostars.org/p/66538/#66569

use strict;
use warnings;
use Bio::Seq;
use Bio::AlignIO;
use Bio::Tools::IUPAC;

my $alnio = Bio::AlignIO->new(-fh => \*DATA, -format => 'fasta');

while (my $aln = $alnio->next_aln) {
    my $cs = $aln->consensus_iupac;
    my $seq = Bio::Seq->new(-seq => $cs, -alphabet => 'dna');
    print "With IUPAC codes: $cs\n";
    print "IUPAC codes meaning: ";
    my $iupac  = Bio::Tools::IUPAC->new(-seq => $seq);
    my %dnasymbols = $iupac->iupac_iub;
    my @residues = split //, $cs;
    for my $r (@residues) {
	if (exists $dnasymbols{$r}) {
	    if (scalar @{$dnasymbols{$r}} > 1) {
		print "[";
		print join "|", @{$dnasymbols{$r}};
		print "]";
	    }
	    else {
		print @{$dnasymbols{$r}};
	    }
	}
    }
    print "\n";
}

__DATA__
>seq1
ACGGGTA
>seq2
GCGGGTC
>seq3
ACGGGTC
>seq4
GCGGGTA
