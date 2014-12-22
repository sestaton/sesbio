#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;
use Getopt::Long;
use File::Basename;

my $infile;
my $outfile;
my $length;
my $percentID;

GetOptions(
	   'i|infile=s'         => \$infile,                    
           'o|outfile=s'        => \$outfile,
	   'l|length=i'         => \$length,
	   'p|percentid=i'      => \$percentID,
	  );
 
if (!$infile && !$outfile) {
    say "\nERROR: No input was given at the command line.\n";
    usage();
    exit(1);
}

my @matches;
my $total    = 0;
$length    //= 0;
$percentID //= 80;

open my $in, '<', $infile or die "\nERROR: Can't open file: $infile\n";
open my $out, '>', $outfile or die "\nERROR: Can't open file: $outfile\n";

while (<$in>) {
    chomp;
    my @hits = split /\t/, $_;
    if ( ($hits[3] >= $length) && ($hits[2] >= $percentID) ) {
	$total += $hits[3];
    }
}
say $out join "\t", $infile, $total;
close $in;
close $out;

exit;
#
# methods
#
sub usage {
    my $script = basename($0);
    print STDERR <<END
USAGE: $script -i inreport -o outreport <options>

Required:
 -i|infile       :     The parsed report the search.
 -o|outfile      :     The summarized output.

Options:
 -l|length      :     Length of match to accept. (Default: no cutoff).
 -p|percentid   :     Percent identity threshold. (Default: 80).\n

END
}
