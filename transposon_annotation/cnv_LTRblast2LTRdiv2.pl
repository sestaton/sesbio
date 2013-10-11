#!/usr/bin/perl -w
#-------------------------------------------------------+
# 
#-------------------------------------------------------+
# DESCRIPTION: 
#                                                       |
# AUTHOR: S. Evan Staton                                |
# CONTACT: statonse<at>uga.edu                          |
# STARTED: 12/15/10
# UPDATED:
#                                                       |
# USAGE: 
#                         
#                                                       |
#-------------------------------------------------------+
# TODO: 
#
# 

use strict;
use Getopt::Long;

my $infile;
my $outfile;
my $usage = "USAGE: cnv_LTRblast2LTRdiv.pl -i inreport -o parsedreport

\tFirst, run the script cng_header2BACname.pl on the file of 5prime sequences
\tand 3prime sequences separately, including the name 5prime and 3prime, respectively. 
\tThis will make it possible to distinguish the five-prime and three-prime LTR sequences 
\tin the blast report (and ignore other hits). 

\tNext, BLAST the LTR sequences together with a low e-value (e.g. 1e-20) and include
\tthe -m 8 option when running blastall. 

\tLast, simply run this script on the BLAST report.";

#counters
#my $total_hits = 0;
#my $parsed_hits = 0;

GetOptions(
           "i|infile=s"           => \$infile,
           "o|outfile=s"          => \$outfile,
          );

# open the infile or die with a usage statement
if ($infile && $outfile) {
    open (INFILE, "<$infile") || print "Error: can't open $infile\n";

    open (OUTFILE, ">$outfile");
}
else {
    if (!$infile){
        die "\n","ERROR: No infile was given at the command line\n\n",$usage,"\n\n"; 
    }
    if (!$outfile){
        die "\n","ERROR: No outfile was given at the command line\n\n",$usage,"\n\n";
    }
}

my @header =  "#5primeLTR"."\t".
              "3primeLTR"."\t".
              "Percent_identity"."\t".
              "Alignment_length"."\n";

print OUTFILE @header;

#my $match_count;

while (<INFILE>) { 
    chomp; 
    my @blfields = split(/\t/, $_);

    my $fiveprime = $blfields[0];
    my $threeprime = $blfields[1];

    if ($blfields[0] !~ $blfields[1]) {

	$fiveprime =~ s/5prime$//;
	$threeprime =~ s/3prime$//;
	  
	#my %seen;

	#if ( ! $seen{$blfields[0].$blfields[1]} ++ ) {

        if ($fiveprime =~ m/$threeprime/) { 
			    	    	 
	    #$match_count++;
	    
	    my @LTRdiv = $blfields[0]."\t".
		         $blfields[1]."\t".
	                 $blfields[2]."\t".
	                 $blfields[3]."\n";

	    print OUTFILE @LTRdiv;
	   
	    #foreach (@LTRdiv) {

	    #	chomp;
	    #	my @out = split(/\t/, $_);
	    # 	my $five = $out[0];
	    #	$five =~ s/5prime$//;
	    #	my $three = $out[1];
	    #	$three =~ s/3prime$//;

       	    #	if ( !$seen{$five.$three}++ ) {
	   	   		    
	    #	    print OUTFILE $_."\n";
		
	    #   }
	    #}

	 }
	
    }

}

exit;
