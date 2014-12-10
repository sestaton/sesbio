#!/usr/bin/env perl

## This script finds all Fastq files in the working directory and
## and writes a library file for use with GapFiller. Currently,
## the match and negation patterns are hardcoded, but this could easily
## be made an option.

use 5.010;
use strict;
use warnings;
use Cwd;
use File::Find;
use List::MoreUtils qw(natatime);
use Getopt::Long;

my $libname;
my $tool;
my $insert_size;
my $deviation;
my $orientation;

GetOptions(
	   'l|libname:s'     => \$libname,
	   't|tool:s'        => \$tool,
	   'd|deviation:f'   => \$deviation,
	   'o|orientation:s' => \$orientation,
	   );

$libname     //= 'lib200';
$tool        //= 'bwa';
$insert_size //= '200';
$deviation   //= '0.25';
$orientation //= 'FR';

my $cwd = getcwd();
my @files;
# there is one negated file because it is unpaired, which won't work with GapFiller
find( sub { 
    push @files, $_ if -f and /^HA0001/ && /\.fq$/ && /[^HA0001_620J2AAXX_1_1_trimmed101.fq]/;
      }, $cwd);
 
my @sorted = sort @files;

my $it = natatime 2, @sorted;
while (my @vals = $it->()) {
    say join q{ }, $libname, $tool, @vals, $insert_size, $deviation, $orientation;
}
