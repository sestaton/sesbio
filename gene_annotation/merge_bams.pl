#!/usr/bin/env perl

use strict;
use warnings;
use File::Find;

my $dir = shift;
my $out = shift;

my @files;
find( sub { push @files, $File::Find::name if -f and /\.bam$/ }, $dir);

if (@files > 0) {
    my $inc_str;
    for my $f (@files) {
	$inc_str .= "I=$f ";
    }
    
    my $cmd = "java -jar ~/apps/picard-tools-2.1.0/picard.jar MergeSamFiles ".
	"CREATE_INDEX=true ".
	"$inc_str ".
	"O=$out";
    system($cmd);
}
