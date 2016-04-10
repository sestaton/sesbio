#!/usr/bin/env perl

## NB: Requires Picard (https://broadinstitute.github.io/picard) to be installed (TODO: Check if this is available)

use strict;
use warnings;
use File::Find;

my $usage = "$0 directory_of_bams output_merged.bam";
my $dir = shift or die $usage;
my $out = shift or die $usage;

my @files;
find( sub { push @files, $File::Find::name if -f and /\.bam$/ }, $dir);

if (@files > 0) {
    my $inc_str;
    for my $f (@files) {
	$inc_str .= "I=$f ";
    }
    
    my $cmd = "java -jar picard.jar MergeSamFiles ".
	"CREATE_INDEX=true ".
	"$inc_str ".
	"O=$out";
    system($cmd);
}
