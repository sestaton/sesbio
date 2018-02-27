#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;
use Bio::AlignIO;

my $usage = "\nUSAGE: get_aln_stats.pl infile\n\n";
my $infile = shift or die $usage;

my $pos = 0;

my $gap_arr = parse_aln($infile);

# NB: per column alignment results are not very useful, better to not print or give
# HSP stats
#for my $hash (@$gap_arr) {
    #for my $key (%$hash) {
	#$pos++;
	#if (defined $hash->{$key}) {
	    #say "Looking at pos: $pos"," which contains: ", $hash->{$key};
	#}   
    #}
#}

exit;
#
# methods
#
sub parse_aln {
    my ($aln_file) = @_;
    my $aln_in = Bio::AlignIO->new(-file   => $aln_file, -format => 'fasta');
				   
    while ( my $aln = $aln_in->next_aln() ) {
	my $percentID = sprintf("%.2f", $aln->percentage_identity);
	my $sim = $percentID/100;
	my $dive = 1-$sim;

	say "align length is: ",$aln->length;
	say "Num residues is: ",$aln->num_residues;
	say "Is flush (bool): ",$aln->is_flush;
	say "Num sequences:   ",$aln->num_sequences;
	say "\nPercent ident:   ",sprintf("%.2f", $aln->percentage_identity);
	say "Similarity:      ",$sim;
	say "Divergence:      ",$dive,"\n";
	#say "Cons string(50): ",$aln->consensus_string(50);
	#say "Cons string: ",$aln->consensus_string;
	return $aln->gap_col_matrix;
    }
}

