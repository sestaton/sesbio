#!/usr/bin/env perl

=head1 NAME 
                                                                       
fasta_to_fastq.pl - Extract the seq and qual string from a Fastq file given a Fasta

=head1 SYNOPSIS    

 fasta_to_fastq.pl -fa seqs.fasta -fq orig_seqs.fastq.gz -o seqs.fastq

=head1 DESCRIPTION
                                                                   
After quality trimming or screening reads one often needs to go back to Fastq for mapping
or assembly. This script takes a Fasta file and will extract the exact sequence and correct
qual string from a Fastq file. The DNA sequences are expected to be different, however, the
sequence names must be the same.

=head1 AUTHOR 

S. Evan Staton                                                

=head1 CONTACT
 
statonse at gmail dot com

=head1 REQUIRED ARGUMENTS

=over 2

=item -fa, --fasta

The Fasta file from which to select sequences.

=item -fq, --fastq

The Fastq file from which to pull the sequence and quality scores.

=item -o, --outfile

The file to write the Fastq records that are selected.

=back

=head1 OPTIONS

=over 2

=item -h, --help

Print a usage statement. 

=item -m, --man

Print the full documentation.

=cut   

use 5.010;
use strict;
use warnings;
use File::Basename;
use autodie qw(open);
use Getopt::Long;
use Pod::Usage

my $fasta;
my $fastq;
my $outfile;
my $help;
my $man;

GetOptions(
	   'fa|fasta=s'   => \$fasta,
	   'fq|fastq=s'   => \$fastq,
	   'o|outfile=s'  => \$outfile,
	   'h|help'       => \$help,
	   'm|man'        => \$man,
	   ) or pod2usage( "Try 'basename($0) --man' for more information." );;

pod2usage( -verbose => 2 ) if $man;
usage() and exit(0) if $help;

if (!$fasta || !$fastq || !$outfile) {
    say "\nERROR: No input was given.";
    usage();
    exit(1);
}

unless (-e $fasta && -e $fastq) {
    say "\nERROR: One or more input files not found.";
    usage();
    exit(1);
}

my $fqct = 0;
my ($fa_idx, $fact) = make_fasta_index($fasta);
my $fq = get_fh($fastq);
open my $out, '>', $outfile;

my ($name, $comm, $seq, $qual);
my @aux = undef;
while (($name, $comm, $seq, $qual) = readfq(\*$fq, \@aux)) {
    if (exists $fa_idx->{$name}) {
	$fqct++;
	my $seq_match;
	if ($seq =~ /($fa_idx->{$name})/) {
	    my ($seq_match, $start, $end) = ($1, $-[0], $+[0]);
	    my $qual_region = substr $qual, $start, length($seq_match);
	    say $out join "\n", "@".$name, $seq_match, q{+}, $qual_region;
	}
    }
}
close $fq;
close $out;

say "$$fact sequences from $fasta were indexed.";
say "$fqct sequences were matched in $fastq and written to $fastq.";

#
# methods
#
sub make_fasta_index {
    my ($fasta) = @_;

    my $fact = 0;
    my $fa = get_fh($fasta);
    my %index;
    my ($name, $comm, $seq, $qual);
    my @aux = undef;
    while (($name, $comm, $seq, $qual) = readfq(\*$fa, \@aux)) {
	$fact++;
	$index{$name} = $seq;
    }
    close $fa;
    return \%index, \$fact;
}

sub get_fh {
    my ($file) = @_;

    my $fh;
    if ($file =~ /\.gz$/) {
        open $fh, '-|', 'zcat', $file or die "\nERROR: Could not open file: $file\n";
    }
    elsif ($file =~ /\.bz2$/) {
        open $fh, '-|', 'bzcat', $file or die "\nERROR: Could not open file: $file\n";
    }
    else {
        open $fh, '<', $file or die "\nERROR: Could not open file: $file\n";
    }

    return $fh;
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
  print STDERR <<END
USAGE: $script -fa seqs.fa -fq seqs.fq -o myseqs.fq 

Required:
    -fa|fasta   :    Fasta file of sequences to map IDs and lengths
    -fq|fastq   :    Fastq file to pull reads from
    -o|outfile  :    The file to place the selected reads.

Options:
    -h|help     :    Print usage statement.
    -m|man      :    Print full documentation.
END
}
