#!/usr/bin/perl -w




use strict;
use File::Basename;
use File::Copy;

my $usage = "\nhelsearch_prepare.pl indir\n";
my $indir = $ARGV[0] || die "\nERROR: No input directory found!\n",$usage;

opendir(DIR,$indir) || die "\nERROR: Could not open directory: $indir\n";
my @fastas = grep /\.fasta$/, readdir DIR;
closedir(DIR);

chdir($indir);

foreach my $fas (@fastas) {
    my ($file,$dir,$ext) = fileparse($fas, qr/\.[^.]*/);
    move("$fas","$file") || die "Copy failed: $!";
}

