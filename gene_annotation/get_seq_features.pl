#!/usr/bin/env perl

use strict;
use warnings;
use Bio::SeqIO;
use Getopt::Long;

my $usage = "perl $0 -i dir -o out\n";
my $indir;
my $outfile;

GetOptions(
           'i|indir=s'   => \$indir,
           'o|outfile=s' => \$outfile,
          );

die $usage if !$indir;

$indir =~ s/\/$// if $indir =~ /\/$/;
my @embl_files = glob("$indir/*.embl");
if (scalar @embl_files < 1) {
    print "\nERROR: Could not find any embl files in $indir. Must end with \".embl\". Exiting.\n";
    exit(1);
}

for my $file (@embl_files) {
    my $seqio = Bio::SeqIO->new(-file => $file, -format => 'EMBL');
    my $seq = $seqio->next_seq;
    for my $feat_obj ($seq->get_SeqFeatures) {
	if ($feat_obj->primary_tag eq "CDS") {
	    print $feat_obj->spliced_seq->seq,"\n";
	    if ($feat_obj->has_tag('gene')) {
		for my $gene ($feat_obj->get_tag_values('gene')) {
		    print "gene: ", $gene,"\n";
		}
	    }
	}
    }
}
