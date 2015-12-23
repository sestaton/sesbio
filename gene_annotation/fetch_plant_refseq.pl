#!/usr/bin/env perl

# consider do this asynchronously with AnyEvent::FTP::Client,
# or in parallel with a (thread) queue

use 5.010;
use strict;
use warnings;
use autodie qw(open);
use Net::FTP;

my $host = "ftp.ncbi.nlm.nih.gov";
my $dir  = "/refseq/release/plant";
my $faa  = shift;
$faa     //= "plant_refseq_all.faa.gz";
open my $outfh, '>>', $faa;

my $ftp = Net::FTP->new($host, Passive => 1, Debug => 0)
    or die "Cannot connect to $host: $@";

$ftp->login or die "Cannot login ", $ftp->message;

$ftp->cwd($dir)
    or die "Cannot change working directory ", $ftp->message;

my @faafiles = grep /\.faa.gz$/, $ftp->ls();
my @sorted   = map  { $_->[0] }
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

    say "=====> working on file: $file";

    collate($file, $outfh);
    unlink $file;
    sleep 1;
}

$ftp->quit;

sub collate {
    my ($file_in, $fh_out) = @_;
    my $lines = do { 
	local $/ = undef; 
	open my $fh_in, '<', $file_in or die "\nERROR: Could not open file: $file_in\n";
	<$fh_in>;
    };
    print $fh_out $lines;
}
