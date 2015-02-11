#!/usr/bin/env perl

## This script finds all Fastq files in the working directory (default), or
## a specified directory and writes a library file for use with GapFiller. 

use 5.010;
use strict;
use warnings;
use Cwd;
use File::Find;
use File::Basename;
use List::MoreUtils qw(natatime);
use Getopt::Long;

my $libname;
my $tool;
my $insert_size;
my $deviation;
my $orientation;
my $directory;
my $match;
my $negate;
my $help;

GetOptions(
	   'l|libname:s'     => \$libname,
	   't|tool:s'        => \$tool,
           'i|insertsize:i'  => \$insert_size,
	   'd|deviation:f'   => \$deviation,
	   'o|orientation:s' => \$orientation,
	   'd|directory:s'   => \$directory,
           'm|match:s'       => \$match,
           'n|negate:s'      => \$negate,
           'h|help'          => \$help,
	   );

usage() and exit(0) if $help;
usage() and exit(1) if !$insert_size || !$libname;

$tool        //= 'bwa';
$deviation   //= '0.25';
$orientation //= 'FR';
$match       //= 'HA0001';
$negate      //= 'HA0001_620J2AAXX_1_1_trimmed101.fq';

my $cwd = getcwd();
$directory //= $cwd
my @files;

# if any file is unpaired it should be negated because these won't work with GapFiller
find( sub { 
    push @files, $_ if -f and /^$match/ && /\.fq$/ && /[^$negate]/;
      }, $directory);
 
my @sorted = sort @files;

my $it = natatime 2, @sorted;
while (my @vals = $it->()) {
    say join q{ }, $libname, $tool, @vals, $insert_size, $deviation, $orientation;
}

sub usage {
    my $script = basename($0);
    print STDERR <<END;

USAGE: $script -l lib200 -i 200

Required:
    -l|libname        :    A library name to be used for the gap filling run.
    -i|insertsize     :    The insert size of the reads.

Options:
    -d|directory      :   The path to search for sequence files.
    -t|tool           :   The name of the toolto use for alignment (Default: bwa).
    -d|deviation      :   The allowed amount of variation between the reads (Default: 0.25).
    -o|orientation    :   The orientation of the reads(Default: FR).
    -h|help           :    Print usage statement.

END
}
