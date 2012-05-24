#!/usr/bin/perl -w

# TODO: 
use strict;
use Bio::AlignIO;
#use Data::Dumper;

my $usage = "\nUSAGE: get_aln_stats.pl infile\n\n";
my $infile = shift or die $usage;

my $pos = 0;

my @gap_arr = parse_aln($infile);

foreach my $hash (@gap_arr) {
    #print Dumper $hash;
    for my $key (%$hash) {
	$pos++;
	if (defined $hash->{$key}) {
	    print "Looking at pos: $pos"," which contains: $hash->{$key}\n";
	}   
    }
}

exit;
#
# Subs
#
sub parse_aln {

    my $aln_file = shift;
    my $aln_in = Bio::AlignIO->new(-file   => $aln_file,
				   -format => 'fasta');
				   
    while ( my $aln = $aln_in->next_aln() ) {

	my $percentID = sprintf("%.2f",$aln->percentage_identity);
	my $sim = $percentID/100;
	my $dive = 1-$sim;

	print "align length is: ",$aln->length, "\n";
	print "Num residues is: ",$aln->no_residues, "\n";
	print "Is flush (bool): ",$aln->is_flush, "\n";
	print "Num sequences:   ",$aln->no_sequences, "\n";
	print "\nPercent ident:   ",sprintf("%.2f",$aln->percentage_identity), "\n";
	print "Similarity:      ",$sim,"\n";
	print "Divergence:      ",$dive,"\n\n";
	#print "Cons string(50): ",$aln->consensus_string(50), "\n";
	#print "Cons string: ",$aln->consensus_string, "\n";
	#print Dumper $aln->gap_col_matrix;
	return $aln->gap_col_matrix;

    }

}

