#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use autodie qw(open);
use lib qw(/home/jmblab/statonse/apps/perlmod/Data-Dump-1.21/blib/lib);
use Data::Dump qw(dd);
use File::Path qw(make_path);
use File::Spec;
use File::Basename;
use POSIX qw(strftime);
use feature 'say';

my $usage = "$0 cls_file fasta_file\n";
my $cls_file = shift or die $usage;
my $fas_file = shift or die $usage;
my $cluster_size;

GetOptions(
	   'c|cls_file=s'  => \$cls_file,
	   'f|fas_file=s'  => \$fas_file,
	   );

open(my $cls, '<', $cls_file);

my $seqhash = fas2hash($fas_file);
#dd $seqhash;

{
    local $/ = '>';
   
    my $cls_dir_path;

    $cluster_size = defined($cluster_size) ? $cluster_size : '500';

    open(my $cls, '<', $cls_file);
    local $/ = '>';

    my $str = POSIX::strftime("%m_%d_%Y_%H_%M_%S", localtime);
    my ($iname, $ipath, $isuffix) = fileparse($cls_file, qr/\.[^.]*/);
    my $cls_dir_base = $iname;
    $cls_dir_base =~ s/\.[^.]+$//;
    my $cls_dir = $cls_dir_base."_cls_fasta_files_$str";
    $cls_dir_path = $ipath.$cls_dir;
    #say "str $str iname $iname cls_dir_base $cls_dir_base cls_dir_path $cls_dir_path";
    make_path($cls_dir_path, {verbose => 0, mode => 0711,});                                                                                                 
    #mkdir($cls_dir_path, 0777) || die "ERROR: Could not make directory for clusters. Exiting.\n";

    while (my $line = <$cls>) {
        my ($clsid, $seqids) = split /\n/, $line;
        next unless defined $seqids;
        my @ids = split /\s+/, $seqids;
        my ($clid, $cls_seq_num) = $clsid =~ /(\S+)\s(\d+)/;
        #$cls_seq_num =~ s/^.*?(\d+)/$1/s;                                                                                                                    
        #say "clid is $clid and cls_seq_num is $cls_seq_num";                                                                                                 
        if ($cls_seq_num >= $cluster_size) {
            say $cls_seq_num;
            my $indiv_clsfile = join "_", $clid, $cls_seq_num;
            $indiv_clsfile .= ".fas";
            #say $indiv_clsfile;                                                                                                                              
            my $cls_file_path = File::Spec->catfile($cls_dir_path, $indiv_clsfile);
            open(my $clsout, '>', $cls_file_path);
            for my $seq (@ids) {
                if (exists $seqhash->{$seq}) {
                    say $clsout join "\n", ">".$seq, $seqhash->{$seq};
                }
            }
            close($clsout);
        }
    }
}

#
# Subs
#

sub fas2hash {
    my $fas_file = shift;
    open(my $fas, '<', $fas_file);

    local $/ = '>';

    my %seqhash;
    while (my $line = <$fas>) {
	my ($seqid, $seq) = split /\n/, $line;;
	$seqhash{$seqid} = $seq;
    }
    close($fas);
    
    return \%seqhash;
}
