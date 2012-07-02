#!/usr/bin/perl -w

use strict;
use Bio::SeqIO;
use Getopt::Long;
use File::Basename;

#
# lexical vars with scope
#
my $infile;
my $outfile;
my $name;
my $format;
my $end;
my $start;
my $desc;
my $help;

GetOptions (
	    'i|infile=s'      => \$infile,
	    'o|outfile=s'     => \$outfile,
	    'n|name=s'        => \$name,
            'sf|seq_format=s' => \$format,
            'end'             => \$end,
	    'start'           => \$start,
	    'description'     => \$desc,
	    'h|help'          => \$help,
	    );

#
# check @ARGV
#
usage() and exit(0) if $help;

if (!$infile || !$outfile || 
    !$name || 
    !$start && !$end) {
    print "\nERROR: No input was given.\n";
    usage();
    exit(1);
}

$format = defined($format) ? $format : 'fasta';

#
# create SeqIO objects
#
my $seq_in  = Bio::SeqIO->new(-format => $format, 
			      -file => $infile); 

my $seq_out = Bio::SeqIO->new(-format => $format,
                              -file => ">$outfile"); 

while( my $seq = $seq_in->next_seq() ) {
    if ($start) {
	my $seqID = $name."_".$seq->id; 
	$seq->id($seqID);
	$seq->description("") unless $desc;
	$seq_out->write_seq($seq);
    }
    if ($end) {
	my $seqID = $seq->id."_".$name;
	$seq->id($seqID);
	$seq->description("") unless $desc;
	$seq_out->write_seq($seq);
    }
}

exit;

#
# subs
#
sub usage {
    my $script = basename($0);
    print STDERR <<END

USAGE: $script -i seqs.fas -o seqs_renamed.fasta -n name [-sf] [--start] [--end] [--description]

Required:
    -i|infile        :    Fasta file to reformat (contig or chromosome).
    -o|outfile       :    File name to write renamed sequences.
    -n|name          :    The name to append to the beginning or end of the primary ID.

Options:
    -sf|seq_format   :    The format of the sequence (Default: fasta).
    --start          :    Place the argument to option [--name] at the beginning of the primary ID.
    --end            :    Place the argument to option [--name] at the end of the primary ID.
    --description    :    If given, the sequence description (anything after the first space in the ID) will be retained
                          (Default: remove anything after the first space).

END
}

