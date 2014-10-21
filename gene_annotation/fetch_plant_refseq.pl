#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;
use autodie qw(open);
use Net::FTP;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);

my $host = "ftp.ncbi.nlm.nih.gov";
my $dir  = "/refseq/release/plant";
my $faa  = "plant_refseq_all.faa";
open my $outfh, '>>', $faa;

my $ftp = Net::FTP->new($host, Passive => 1, Debug => 0)
    or die "Cannot connect to $host: $@";

$ftp->login or die "Cannot login ", $ftp->message;

$ftp->cwd($dir)
    or die "Cannot change working directory ", $ftp->message;

my @faafiles = grep /\.faa.gz$/, $ftp->ls();
my @sorted = map  { $_->[0] }
             sort { $a->[1] <=> $b->[1] }
             map  { [ $_, /(\d+)/ ] }
             @faafiles;

for my $file (@sorted) {
    $ftp->binary();
    my $rsize = $ftp->size($file) or die "Could not get size ", $ftp->message;
    $ftp->get($file) or die "get failed ", $ftp->message;
    my $lsize = -s $file;
    die "Failed to fetch complete file: $file (local size: $lsize, remote size: $rsize)"
	unless $rsize == $lsize;
    my $flatdb = $file;
    $flatdb =~ s/\.gz$//;

    say "=====> working on file: $file";
    my $status = gunzip $file => $flatdb
        or die "gunzip failed: $GunzipError\n";

    collate($flatdb, $outfh);
    unlink $file;
    unlink $flatdb;
    sleep 1;
}

$ftp->quit;

sub collate {
    my ($flatdb, $outfh) = @_;

    open my $infh, '<', $flatdb;
    while (<$infh>) {
	print $outfh $_;
    }
    close $infh;
 }
    
