#!/usr/bin/perl -w 
#_______________________________________________________________________+
#                                                                       |
# count_bt_bases.pl               
#_______________________________________________________________________+
#                                                                       |
# Description: 
#                                                                       |
# Author: S. Evan Staton                                                |
# Contact: statonse<at>uga.edu                                          |
# Started: 2.22.11                                                      |                                                
# Updated:                                                      
#                                                                       |
# Suggested usage:                                                      |

#_______________________________________________________________________+
# TODO: 

use strict;
use Getopt::Long;
use File::Basename;
#use Statistics::Descriptive;  # just getting sum right now and we don't need help for that ;)
#use File::Copy;               # for copying sequences to their final destination
#use File::Spec;

my $script = basename($0);
# define vars
my $usage = "USAGE: $script -i inreport -o outreport <options>

\tRequired:
\t -i|infile       :     The parsed report the search.
\t -o|outfile      :     The summarized output. 

\tOptions:
\t -l|length      :     Length of match to accept. (Default: no cutoff).
\t -p|percentid   :     Percent identity threshold. (Default: 80).\n";

my $infile;
my $outfile;
my $length;
my $percentID;

GetOptions(# Required arguments
	   "i|infile=s"         => \$infile,                    
           "o|outfile=s"        => \$outfile,                   # make this the alignment report       

	   "l|length=i"         => \$length,
	   "p|percentid=i"      => \$percentID,
	  );
# 
if (!$infile && !$outfile) {
    die "\nERROR: No input was given at the command line.\n\n",$usage,"\n\n";
}

my $length_threshold = defined($length) ? $length : "0";
my $percentID_threshold = defined($percentID) ? $percentID : "80";

my ($file, $dir, $ext) = fileparse($infile, qr/\.[^.]*/);
open( my $in, '<' $infile ) or die "\nERROR: Can't open file: $infile\n";
open( my $out, '>', $outfile ) or die "\nERROR: Can't open file: $outfile\n";

my @matches;
my $total = 0;

while (<$in>) {
    chomp;
    next if /^#/;
    my @hits = split(/\t/,$_);
    if ( ($hits[3] >= $length_threshold) && ($hits[2] >= $percentID_threshold) ) {
	push(@matches,$hits[3]);
	foreach my $match (@matches) {
	    $total += $match;
	}
	print $_."\n";
    }
}
print $out join("\t",($file,$total)),"\n";
close($in);
close($out);

exit;
