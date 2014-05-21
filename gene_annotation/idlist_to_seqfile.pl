#!/usr/bin/env perl

## NB: faSomeRecords from Kent source is the fastest for this task.

use 5.010;
use strict;
use warnings;
use File::Basename;
use Getopt::Long;

# lexical vars
my $idlist;
my $infile;
my $outfile;
my $format;

GetOptions(
	   'id|idlist=s' => \$idlist,
	   'i|infile=s'  => \$infile,
	   'o|outfile=s' => \$outfile,
           'f|format=s'  => \$format,
	   );

# check @ARGV
if (!$infile || !$outfile || !$idlist) {
    usage();
    exit(1);
}

$format //= 'fasta';

if ($format !~ /fasta/i || $format !~ /fastq/i) {
    say "\nERROR: $format not recognized. Must be 'fasta' or 'fastq'. Exiting.\n";
    exit(1);
}

open my $out, '>', $outfile or die "\nERROR: Could not open file: $outfile\n";

my $idlist_hash = id2hash($idlist);
my $seqhash = seq2hash($infile, $format);

for my $gene (keys %$idlist_hash) {
    if (exists $seqhash->{$gene}) {
	if ($format =~ /fasta/i) {
	    say $out join "\n", ">".$gene, $fa_hash->{$gene};
	}
	elsif ($format =~ /fastq/i) {
	    my ($seq, $qual) = split '~~', $seqhash->{$gene};
	    say $out join "\n", "@".$gene, $seq, "+", $qual;
	}
    }
}
close $out;

exit;
#
# methods
#
sub id2hash {
    my $idlist = shift;
    open my $fh, '<', $idlist or die "\nERROR: Could not open file: $!\n";

    my %hash;
    while (<$fh>) {
	chomp;
	$hash{$_} = 1;
    }
    close $fh;
    return \%hash;
}

sub seq2hash {
    my ($infile, $format) = @_;
    open my $fh, '<', $infile or die "\nERROR: Could not open file: $infile\n";

    my @aux = undef;
    my ($name, $comm, $seq, $qual);
    my %seqhash;
    my $seqct = 0;
    
    while (($name, $comm, $seq, $qual) = readfq(\*$fh, \@aux)) {
	$seqct++;
	if ($format =~ /fasta/i) {
	    $seqhash{$name} = $seq;
	}
	elsif ($format =~ /fastq/i) {
	    $seqhash{$name} = join "~~", $seq, $qual;
	}
    }
    
    close $fh;
    return \%seqhash, \$seqct;
}

sub readfq {
    my ($fh, $aux) = @_;
    @$aux = [undef, 0] if (!@$aux);
    return if ($aux->[1]);
    if (!defined($aux->[0])) {
        while (<$fh>) {
            chomp;
            if (substr($_, 0, 1) eq '>' || substr($_, 0, 1) eq '@') {
                $aux->[0] = $_;
                last;
            }
        }
        if (!defined($aux->[0])) {
            $aux->[1] = 1;
            return;
        }
    }
    my ($name, $comm);
    defined $_ && do {
        ($name, $comm) = /^.(\S+)(?:\s+)(\S+)/ ? ($1, $2) : 
	                 /^.(\S+)/ ? ($1, '') : ('', '');
    };
    my $seq = '';
    my $c;
    $aux->[0] = undef;
    while (<$fh>) {
        chomp;
        $c = substr($_, 0, 1);
        last if ($c eq '>' || $c eq '@' || $c eq '+');
        $seq .= $_;
    }
    $aux->[0] = $_;
    $aux->[1] = 1 if (!defined($aux->[0]));
    return ($name, $comm, $seq) if ($c ne '+');
    my $qual = '';
    while (<$fh>) {
        chomp;
        $qual .= $_;
        if (length($qual) >= length($seq)) {
            $aux->[0] = undef;
            return ($name, $comm, $seq, $qual);
        }
    }
    $aux->[1] = 1;
    return ($name, $seq);
}

sub usage {
    my $script = basename($0);
    print STDERR<<EOF
USAGE: $script [-id] [-i] [-o] [-f]

Required:
    -id|idlist    :      An ID list of records (one per line).  
    -i|infile     :      A fasta/q file to pull sequences from.
    -o|outfile    :      A file to write the records to.
    -f|format     :      Input format, must be 'fasta' or 'fastq' (Default: fasta).

EOF
}
