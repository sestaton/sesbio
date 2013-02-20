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
use feature 'say';

my $usage = "$0 cls_file fasta_file\n";
my $cls_file = shift or die $usage;
my $fas_file = shift or die $usage;

GetOptions(
	   'c|cls_file=s'  => \$cls_file,
	   'f|fas_file=s'  => \$fas_file,
	   );

open(my $cls, '<', $cls_file);

my $seqhash = fas2hash($fas_file);
#dd $seqhash;

{
    local $/ = '>';
   
    while (my $line = <$cls>) {
	my ($clsid, $seqids) = split /\n/, $line;
	next unless defined $seqids;
	my @ids = split /\s+/, $seqids;
	my $clsfile = $clsid;
	$clsfile =~ s/\s/\_/g;
	$clsfile =~ s/^\>//;
	$clsfile .= ".fas";
	my ($iname, $ipath, $isuffix) = fileparse($cls_file, qr/\.[^.]*/);
	my $cls_dir = $iname."_cls_cluster_fasta_files";
	my $cls_dir_path = $ipath.$cls_dir;
	make_path($cls_dir_path, {verbose => 0, mode => 0711,});
	my $cls_file_path = File::Spec->catfile($cls_dir_path, $clsfile);
	#say $cls_file_path;
        open(my $clsout, '>', $cls_file_path);
	for my $seq (@ids) {
	    if (exists $seqhash->{$seq}) {
		#say "$seq $seqhash->{$seq}";
		say $clsout join "\n", ">".$seq, $seqhash->{$seq};
	    }
	}
	close($clsout);
    }
}

#
# Subs
#

sub fas2hash {
    my $fas_file = shift;
    open(my $fas, '<', $fas_file);

    local $/ = '>';

    while (my $line = <$fas>) {
	my ($seqid, $seq) = split /\n/, $line;;
	$seqhash{$seqid} = $seq;
    }
    close($fas);
    
    return \%seqhash;
}
