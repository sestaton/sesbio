#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;
use Bio::AlignIO;
use Getopt::Long;

my $infile;
my $outfile;
my $iformat; ## Optional (Default: fasta)
my $oformat; ## Optional (Default: clustalw)

my $usage = "USAGE: perl $0 -i infile -o outfile -if format -of format\n\nOptions 'if' and 'of' are optional (Defaults: fasta for 'if' and clustalw for 'of')."; 

GetOptions(
           'i|infile=s'     => \$infile,
           'o|outfile=s'    => \$outfile,
           'if|informat:s'  => \$iformat,
           'of|outformat:s' => \$oformat,
    );

say $usage and exit(0) if !$infile or !$outfile;

$iformat //= 'fasta';
$oformat //= 'clustalw';

my $alignin  = Bio::AlignIO->new(-file => $infile, -format => $iformat);
my $alignout = Bio::AlignIO->new(-file => ">$outfile", -format => $oformat);
# note: we quote -format to keep older perl's from complaining.

while ( my $aln = $alignin->next_aln() ) {
    $alignout->write_aln($aln);
}
