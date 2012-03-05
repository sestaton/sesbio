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
# VARIABLE SCOPE 
#
my $infile; 
my $outfile;
my $format;
my $tophit;
my $toponly;
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
	   "toponly"       => \$toponly,
	   "l|length=i"    => \$length_thresh,
	   "p|percentid=f" => \$pid_thresh,
	   "v|verbose"     => \$verbose,
	   );
	 
#
# Check @ARGVs 
#
if (!$infile || !$outfile){
    print "\nERROR: No infile was given at the command line\n\n";
    &usage;
    exit(1);
}

if (!$format){
    print "\nERROR: No blastformat was given at the command line\n\n";
}
if ($verbose && !$tophit) {
    print "\nIgnoring option 'verbose' without option 'top'\n\n";
}

my ($tmp_bl_rep) = length_pid_filter($infile);

if ($tophit) {
    get_top_hit($infile,$tmp_bl_rep,$outfile);
    unlink($tmp_bl_rep) if $toponly;
} else {
    move("$tmp_bl_rep","$outfile");
}

exit;

#
# Subroutines
# 
sub length_pid_filter {

    my $bl = shift;

    my $dir = cwd();
    #$dir = $dir."/";
    my $tmp = File::Temp->new( TEMPLATE => $bl.'tempblastrepXXXX',
    #			       DIR => $dir,
    			       UNLINK => 0);
    my $temp_fname = $tmp->filename;
    #my ($tmp_fh, $temp_fname) = tempfile($bl."tempblastrepXXXX",
    #					 DIR => $dir,
    #					 UNLINK => 0);
					 
    open(my $allout, '>', $temp_fname) or die "\nERROR: Could not open file: $temp_fname\n";

    # create SearchIO object to read in blast report and to write outfile
    my $search_in = Bio::SearchIO->new(-format => $format, -file => $bl);

    my $header = "#Query\tTotal_hits\tHit\tHSP_Length\tPercent_ID\tPercent_Cov\tHit_Description\tHit_Significance\n";
    print $allout $header;

    while ( my $result = $search_in->next_result ) {

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
			    
			    print $allout $bl_hits;
			    
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
			
			print $allout $bl_hits;
			
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
			
			print $allout $bl_hits;

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
		    
		    print $allout $bl_hits;

		}
	    }
	}
    }

    close($allout);
    return($temp_fname);

}

sub get_top_hit {

    my ($bl_in,$hits,$bl_out) = @_;

    open(my $allin , '<', $hits) or die "\nERROR: Could not open file: $hits\n";
    open(my $topout , '>', $bl_out) or die "\nERROR: Could not open file: $bl_out\n";

    my $allq = 0;
    my $uniqq = 0;
    my %seen = ();

    while (<$allin>) { 
	chomp;
	if (/^#/) {
	    print $topout $_."\n";
	} else {
	#unless (/^#/) {
	    $allq++;
	    my @blfields = split(/\t/, $_);
	    my $key = $blfields[0];
	    if (! $seen{ $key }++) {
		print $topout $_."\n";
	    }
	}
    }

    $uniqq += keys %seen;

    if ($verbose) {
	print "$allq Total sequences in report $bl_in.\n";
	print "$uniqq Unique sequences have been written to: $bl_out\n";
    }

    close($allin);
    close($topout);  

}

sub usage {

    my $script = basename($0);
        print STDERR <<END
USAGE: $script [-i][-o][-f][-t][--toponly][-l][-p][-v]

Required:
   -i|infile     :     The BLAST report to parse.
   -f|format     :     The BLAST format (e.g. blastxml, blast);
   -o|outfile    :     A file to place the desired BLAST results.

Options:
   -t|top        :     Print the top BLAST hit for each query sequence.
   --toponly     :     Print only the top BLAST hit for each query sequence,
                       and discard all the other hits.
   -l|length     :     Keep only hits above a certain length threshhold (integer).
   -p|percentid  :     Keep only hits above a certain percent identity threshhold (integer).
   -v|verbose    :     Print information about the number of BLAST hits (ignored unless
                       used in conjunction with --top option).

END
}
