#!/usr/bin/perl -w

# Convert a blastxml report to a blast table
     
# 
#TODO: 
#
use strict;
use Bio::SearchIO; # for searchio stuff with blast reports
use Getopt::Long; # for getting options at the command line
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
if (!$infile){
    print "\nERROR: No infile was given at the command line\n\n";
    usage();
    exit(0);
}
if (!$outfile){
    print "\nERROR: No outfile was given at the command line\n\n";
    usage();
    exit(0);
}
if (!$format){
    print "\nERROR: No blastformat was given at the command line\n\n";
}

my $allq = 0;
$length_thresh = defined($length_thresh) ? $length_thresh : '0';          # we are going to set defaults this way
$pid_thresh = defined($pid_thresh) ? $pid_thresh : '0';                   # to work with Perl versions released prior to 5.10

open(my $blastout, '>', $outfile) or die "\nERROR: Could not open file: $!\n";

# create SearchIO object to read in blast report and to write outfile
my $search_in = Bio::SearchIO->new(-format => $format, -file => $infile, -tempfile => 1);


#my $header = "#Query\tHit\tPercent_ID\tHSP_len\tNum_mismatch\tNum_gaps\tQuery_start\tQuery_end\tHit_start\tHit_end\tE-value\tBit_score\n";
#print $blastout $header;

while ( my $result = $search_in->next_result ) {
    $allq++;
    while( my $hit = $result->next_hit ) {
	my $hitlen = $hit->length();
	    
	while( my $hsp = $hit->next_hsp ) {
	    my $hsplen = $hsp->length('total');
	    
	    if( $hsplen >= $length_thresh && $hsp->percent_identity >= $pid_thresh ) {
			
		my $percent_identity = sprintf("%.2f",$hsp->percent_identity);
		
		my @matches = $hsp->matches('hit');
		my $matches = @matches;
		my $mismatches = $hsplen - $matches;
		
		my $bl_hits = $result->query_name."\t".
		    $hit->name."\t".
		    $percent_identity."\t".
		    $hsplen."\t".
		    $mismatches."\t".
		    $hsp->gaps."\t".
		    $hsp->start('query')."\t".
		    $hsp->end('query')."\t".
		    $hsp->start('hit')."\t".
		    $hsp->end('hit')."\t".
		    $hsp->evalue."\t".
		    $hsp->bits."\n";
       		    			
		print $blastout $bl_hits;
		
	    }  
	}
    }
}

close($blastout);

if ($verbose) {
    print "$allq total hits written to report $blastout.\n";
}

sub usage {

    my $script = basename($0);
        print STDERR <<END
USAGE: $script [-i][-o][-f][-t][-l][-p][-v]

Required:
   -i|infile     :     The BLAST report to parse.
   -f|format     :     The BLAST format (e.g. blastxml, blast);
   -o|outfile    :     A file to place the desired BLAST results.

       Options:
   -t|top        :     Print the top BLAST hit for each query sequence.
   -l|length     :     Keep only hits above a certain length threshhold (integer).
   -p|percentid  :     Keep only hits above a certain percent identity threshhold (integer).
   -v|verbose    :     Print information about the number of BLAST hits (ignored unless
                       used in conjunction with --top option).

NB: The output format is identical to the blasttable format (i.e., "-m 8" with legacy BLAST).

END
}
