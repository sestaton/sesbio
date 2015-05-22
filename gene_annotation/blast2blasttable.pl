#!/usr/bin/env perl

# Convert a blastxml or standard blast report to a blast table
# TODO: Not currently getting top hit.

use 5.010;
use strict;
use warnings;
use Bio::SearchIO;
use Getopt::Long;
use File::Basename;

#
# Vars
#
my $infile; 
my $outfile;
my $format;
my $tophit;
my $length_thresh;
my $pid_thresh;
my $verbose;

#
# OPTIONS  
#
GetOptions(
	   'i|infile=s'    => \$infile,
	   'o|outfile=s'   => \$outfile,
	   'f|format=s'    => \$format,
	   't|top'         => \$tophit,
	   'l|length=i'    => \$length_thresh,
	   'p|percentid=f' => \$pid_thresh,
	   'v|verbose'     => \$verbose,
	   );
	 
#
# Check @ARGVs 
#
if (!$infile || !$outfile || !$format){
    say "\nERROR: Command line not parsed correctly.\n";
    usage();
    exit(0);
}

my $allq = 0;
$length_thresh //= 0;
$pid_thresh //= 0;

open my $blastout, '>', $outfile or die "\nERROR: Could not open file: $!\n";

# create SearchIO object to read in blast report and to write outfile
my $search_in; 
if ($tophit) {
    $search_in = Bio::SearchIO->new(-format => $format, -file => $infile, -tempfile => 1, -best_hit_only => 1);
}
else {
    $search_in = Bio::SearchIO->new(-format => $format, -file => $infile, -tempfile => 1);
}
say STDERR "Getting top hits only..." if $tophit;
#say $blastout join "\t", "#Query", "Hit", "Percent_ID", "HSP_len", "Num_mismatch", "Num_gaps", 
#                          "Query_start", "Query_end", "Hit_start", "Hit_end", "E-value", "Bit_score";

while ( my $result = $search_in->next_result ) {
    while ( my $hit = $result->next_hit ) {
	while ( my $hsp = $hit->next_hsp ) {
	    my $hsplen = $hsp->length('total');
	    if ( $hsplen >= $length_thresh && $hsp->percent_identity >= $pid_thresh ) {
		$allq++;
		my $percent_identity = sprintf("%.2f",$hsp->percent_identity);
		
		my @matches = $hsp->matches('hit');
		my $matches = @matches;
		my $mismatches = $hsplen - $matches;
		
		my $name = $result->query_description;
		$name =~ s/\s+/_/g;
		say $blastout join "\t", $name, $hit->name, 
		$percent_identity, $hsplen, $mismatches, $hsp->gaps,
		$hsp->start('query'), $hsp->end('query'), $hsp->start('hit'), 
		$hsp->end('hit'), $hsp->evalue, $hsp->bits, $hit->description;
	    }  
	}
    }
}
close $blastout;

say "$allq total hits written to report $outfile" if $verbose;

sub usage {
    my $script = basename($0);
        print STDERR <<END
USAGE: $script [-i][-o][-f][-t][-l][-p][-v]

Required:
   -i|infile     :     The BLAST report to parse.
   -f|format     :     The BLAST format (e.g. blastxml, blast);
   -o|outfile    :     A file to place the desired BLAST results.

Options:
   -t|top        :     Print the top BLAST hit for each query sequence (Not Implemented).
   -l|length     :     Keep only hits above a certain length threshhold (integer).
   -p|percentid  :     Keep only hits above a certain percent identity threshhold (integer).
   -v|verbose    :     Print information about the number of BLAST hits..

NB: The output format is identical to the blasttable format 
    (i.e., "-m 8" with legacy BLAST or -outfmt 6 with BLAST+).

END
}
