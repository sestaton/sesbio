#!/usr/bin/perl -w
     
# 
#TODO: 
#
use strict;
use Bio::SearchIO; # for searchio stuff with blast reports
use Getopt::Long; # for getting options at the command line
use File::Copy;
use File::Basename;
use File::Temp qw(tempfile);
use Cwd;

#
# 
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
	   "i|infile=s"    => \$infile,
	   "o|outfile=s"   => \$outfile,
	   "f|format=s"    => \$format,
	   "t|top"         => \$tophit,
	   "l|length=i"    => \$length_thresh,
	   "p|percentid=f" => \$pid_thresh,
	   "v|verbose"     => \$verbose,
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

open(my $blastout, '>', $outfile) or die "\nERROR: Could not open file: $outfile\n";

# create SearchIO object to read in blast report and to write outfile
my $search_in;
if ($tophit) {
    $search_in = Bio::SearchIO->new(-format => $format, -file => $infile, -tempfile => 1, -best => 1);
} else {
    $search_in = Bio::SearchIO->new(-format => $format, -file => $infile, -tempfile => 1);
}

#$search_in->best_hit_only if $tophit;

my $header = "#Query\tTotal_hits\tHit\tHSP_Length\tPercent_ID\tPercent_Cov\tHit_Description\tHit_Significance\n";
print $blastout $header;

while ( my $result = $search_in->next_result ) {
    $allq++;
    while( my $hit = $result->next_hit ) {
	my $hitlen = $hit->length();
	    
	while( my $hsp = $hit->next_hsp ) {
	    my $hsplen = $hsp->length('total');
	    
	    if ($length_thresh && $pid_thresh) {
		
		if( $hsp->length('total') > $length_thresh ) {
		    
		    if ( $hsp->percent_identity >= $pid_thresh ) {
			
			my $percent_identity = sprintf("%.2f",$hsp->percent_identity);
			my $percent_coverage = sprintf("%.2f",$hsplen/$hitlen); 
			
			my $bl_hits = $result->query_name."\t".
			    $result->num_hits."\t".
			    $hit->name."\t".
			    $hsp->length('total')."\t".
			    $percent_identity."\t".
			    $percent_coverage."\t".
			    $hit->description."\t".
			    $hit->significance."\n";
			
			print $blastout $bl_hits;
			
		    }  
		}
	    }
	    if ($length_thresh && !$pid_thresh) {
		
		if( $hsp->length('total') > $length_thresh ) {
		    
		    my $percent_identity = sprintf("%.2f",$hsp->percent_identity);
		    my $percent_coverage = sprintf("%.2f",$hsplen/$hitlen);
		    
		    my $bl_hits = $result->query_name."\t".
			$result->num_hits."\t".
			$hit->name."\t".
			$hsp->length('total')."\t".
			$percent_identity."\t".
			$percent_coverage."\t".
			$hit->description."\t".
			$hit->significance."\n";
		    
		    print $blastout $bl_hits;
		    
		}
	    }
	    if (!$length_thresh && $pid_thresh) {
		
		if ( $hsp->percent_identity >= $pid_thresh ) {
		    
		    my $percent_identity = sprintf("%.2f",$hsp->percent_identity);
		    my $percent_coverage = sprintf("%.2f",$hsplen/$hitlen);
		    
		    my $bl_hits = $result->query_name."\t".
			$result->num_hits."\t".
			$hit->name."\t".
			$hsp->length('total')."\t".
			$percent_identity."\t".
			$percent_coverage."\t".
			$hit->description."\t".
			$hit->significance."\n";
		    
		    print $blastout $bl_hits;
		    
		}
	    } 
	    if (!$length_thresh && !$pid_thresh) {
		
		my $percent_identity = sprintf("%.2f",$hsp->percent_identity);
		my $percent_coverage = sprintf("%.2f",$hsplen/$hitlen);
		
		my $bl_hits = $result->query_name."\t".
		    $result->num_hits."\t".
		    $hit->name."\t".
		    $hsp->length('total')."\t".
		    $percent_identity."\t".
		    $percent_coverage."\t".
		    $hit->description."\t".
		    $hit->significance."\n";
		
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

END
}
