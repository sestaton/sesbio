#!/usr/bin/env perl

## This script takes a cluster file from Transposome and a sequence file (with matching IDs
## as in the cluster file) and generates a FASTA file for each cluster.

use 5.010;
use strict;
use warnings;
use autodie;
use Getopt::Long;
use File::Spec;
use File::Basename;
use File::Path qw(make_path);
use POSIX      qw(strftime);

my $usage = "USAGE: perl ".basename($0)." -c cls_file -f fasta_file [-s cluster_size (integer)]\n";
my $cls_file;
my $fas_file;
my $cluster_size;

GetOptions(
	   'c|cls_file=s'     => \$cls_file,
	   'f|fas_file=s'     => \$fas_file,
           's|cluster_size=i' => \$cluster_size,
	   );

die $usage if !$cls_file || !$fas_file;

my $seqhash = fas2hash($fas_file);
$cluster_size //= 500;

{
    local $/ = '>';
   
    my $cls_dir_path;

    open my $cls, '<', $cls_file;

    my $str = strftime("%m_%d_%Y_%H_%M_%S", localtime);
    my ($iname, $ipath, $isuffix) = fileparse($cls_file, qr/\.[^.]*/);
    my $cls_dir_base = $iname;
    $cls_dir_base =~ s/\.[^.]+$//;
    my $cls_dir = $cls_dir_base."_cls_fasta_files_$str";
    $cls_dir_path = $ipath.$cls_dir;
    make_path($cls_dir_path, {verbose => 0, mode => 0711,});
 
    while (my $line = <$cls>) {
        my ($clsid, $seqids) = split /\n/, $line;
        next unless defined $seqids;
        my @ids = split /\s+/, $seqids;
        my ($clid, $cls_seq_num) = $clsid =~ /(\S+)\s(\d+)/;
        #$cls_seq_num =~ s/^.*?(\d+)/$1/s;                                                                               
        if ($cls_seq_num >= $cluster_size) {
            my $indiv_clsfile = join "_", $clid, $cls_seq_num;
            $indiv_clsfile .= ".fas";
   
            my $cls_file_path = File::Spec->catfile($cls_dir_path, $indiv_clsfile);
            open my $clsout, '>', $cls_file_path;
            for my $seq (@ids) {
                if (exists $seqhash->{$seq}) {
                    say $clsout join "\n", ">".$seq, $seqhash->{$seq};
                }
            }
            close $clsout;
        }
    }
}

#
# methods
#
sub fas2hash {
    my $fas_file = shift;
    open my $fas, '<', $fas_file;

    local $/ = '>';

    my %seqhash;
    while (my $line = <$fas>) {
	my ($seqid, $seq) = split /\n/, $line;
	$seqhash{$seqid} = $seq;
    }
    close $fas;
    
    return \%seqhash;
}
