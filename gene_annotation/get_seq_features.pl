#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;
use Bio::Seq;
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
    say "\nERROR: Could not find any embl files in $indir. Must end with \".embl\". Exiting.";
    exit(1);
}

my $contig;
for my $file (@embl_files) {
    my $seqio = Bio::SeqIO->new(-file => $file, -format => 'EMBL');
    my $seq_obj = $seqio->next_seq;
    my $id = $seq_obj->id;
    $id =~ s/.*annotations\.//;
    $id =~ s/\.final\.$//;
    #my $contig = get_contig_name($seq_obj, $file);
    #die "The file is: $file\n" unless defined $contig;
    for my $feat_obj ($seq_obj->get_SeqFeatures) {
	if ($feat_obj->primary_tag eq "gene") {
	    if ($feat_obj->has_tag('gene')) {
		for my $gene ($feat_obj->get_tag_values('gene')) {
		    my $seq = $feat_obj->spliced_seq->seq;
		    my $seqid;
		    eval { $seqid = $id."_".$gene."_".$feat_obj->location->start."_".$feat_obj->location->end; };
		    if ($@) { say "In $file, could not resolve contig name."; }
		    my $geneseq_file = $seqid;
		    $geneseq_file .= ".fasta";
		    my $geneseq = Bio::Seq->new(-seq => $seq, -id => $seqid);
		    my $geneseqio = Bio::SeqIO->new(-file => ">$geneseq_file", -format => 'fasta');
		    $geneseqio->write_seq($geneseq);
		}
	    }
	}
    }
}


sub get_contig_name {
    my ($seq_obj, $file) = @_;
    my $contig;
    for my $feat_obj ($seq_obj->get_SeqFeatures) {
        if ($feat_obj->primary_tag eq "contig") {
            if ($feat_obj->has_tag('note')) {
                for my $contig_seq ($feat_obj->get_tag_values('note')) {
                    $contig_seq =~ s/.*annotations\.//;
                    $contig_seq =~ s/\.final\.$//;
		    $contig = $contig_seq;
                    unless (defined $contig) {
                        say "problem with $file";
                        exit(1);
                    }
                }
            }
        }
    }
    return $contig;
}
	
