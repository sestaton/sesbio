#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename;
use File::Copy qw(move);

my $usage = "\nhelsearch_prepare.pl indir\n";
my $indir = shift or die "\nERROR: No input directory found!\n",$usage;

opendir my $dir, $indir or die "\nERROR: Could not open directory: $indir\n";
my @fastas = grep /\.fasta$/, readdir $dir;
closedir $dir;

chdir $indir;

for my $fas (@fastas) {
    my ($file, $dir, $ext) = fileparse($fas, qr/\.[^.]*/);
    move($fas, $file) || die "Copy failed: $!";
}

