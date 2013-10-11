#!/usr/bin/env perl
# TODO: 

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
           'o|outfile=s'        => \$outfile,                   # make this the alignment report       
	   'l|length=i'         => \$length,
	   'p|percentid=i'      => \$percentID,
	  );
 
if (!$infile && !$outfile) {
    say "\nERROR: No input was given at the command line.\n";
    usage();
    exit(1);
}

my $length_threshold = defined($length) ? $length : "0";
my $percentID_threshold = defined($percentID) ? $percentID : "80";

#my ($file,$dir,$ext) = fileparse($infile, qr/\.[^.]*/);
open my $in, '<', $infile or die "\nERROR: Can't open file: $infile\n";
open my $out, '>', $outfile or die "\nERROR: Can't open file: $outfile\n";

my @matches;
my $total = 0;

while (<$in>) {
    chomp;
    my @hits = split /\t/, $_;
    if ( ($hits[3] >= $length_threshold) && ($hits[2] >= $percentID_threshold) ) {
	$total += $hits[3];
    }
}
say $out join "\t", $infile, $total;
close $in;
close $out;

exit;
#
# subs
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
