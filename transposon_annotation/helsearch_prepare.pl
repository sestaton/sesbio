#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename;
use File::Copy;

my $usage = "\nhelsearch_prepare.pl indir\n";
my $indir = shift or die "\nERROR: No input directory found!\n",$usage;

opendir(DIR,$indir) || die "\nERROR: Could not open directory: $indir\n";
my @fastas = grep /\.fasta$/, readdir DIR;
closedir(DIR);

chdir($indir);

for my $fas (@fastas) {
    my ($file,$dir,$ext) = fileparse($fas, qr/\.[^.]*/);
    move("$fas","$file") || die "Copy failed: $!";
}

