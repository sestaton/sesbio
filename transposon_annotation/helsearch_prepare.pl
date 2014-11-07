#!/usr/bin/env perl

# This script creates input files in the format required
# by the Helitron-finding program Helsearch.

# Usage: Copy all your Fasta files to a directory, run this script with
#        that directory as the only argument and this script
#        will format the files for input to Helsearch.

use strict;
use warnings;
use File::Basename;
use File::Copy qw(move);

my $usage = "\nhelsearch_prepare.pl indir\n";
my $indir = shift or die "\nERROR: No input directory found!\n",$usage;

opendir my $dir, $indir or die "\nERROR: Could not open directory: $indir\n";
my @fastas = grep /\.fa*/, readdir $dir;
closedir $dir;

chdir $indir;

for my $fas (@fastas) {
    my ($file, $dir, $ext) = fileparse($fas, qr/\.[^.]*/);
    move($fas, $file) || die "Copy failed: $!";
}

