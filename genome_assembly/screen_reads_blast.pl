#!/usr/bin/env perl

=head1 NAME 
                                                                       
screen_reads_blast.pl - Retain portions of query sequences matching a target.

=head1 SYNOPSIS    

screen_reads_blast.pl -f seqs.fas -b seqs_db.bln -o screened_seqs.fas

=head1 DESCRIPTION
                                                                   
Filter matches of a query against a target, and retain only the matching portion over
a certain length.
                                                                       
=head1 AUTHOR 

S. Evan Staton                                                

=head1 CONTACT
 
statonse at gmail dot com

=head1 REQUIRED ARGUMENTS

=over 2

=item -b, --blastfile

The blast file (tab-delimited "-m8" or "-outfmt 6").

=item -i, --infile

The file of reads used as the query in the blast search.

=item -o, --outfile

A file to place the filtered reads.

=back

=head1 OPTIONS

=over 2

=item -l, --length

A minimum length threshold (Default: 50).

=item -h, --help

Print a usage statement. 

=item -m, --man

Print the full documentation.

=cut

use 5.010;
use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use autodie qw(open);

my $blast;
my $fasta;
my $length;
my $outfile;
my $help;
my $man;

GetOptions(
           'i|infile=s'      => \$fasta,
           'b|blastfile=s'   => \$blast,
           'o|outfile=s'     => \$outfile,
           'l|length=i'      => \$length,
           'h|help'          => \$help,
           'm|man'           => \$man,
    ) || pod2usage( "Try 'basename($0) --man' for more information." );

pod2usage( -verbose => 2 ) if $man;
usage() and exit(0) if $help;

if (!$fasta || !$blast || !$outfile) {
    usage();
    exit(0);
}

$length //= 50;
open my $in, '<', $blast;
open my $fas, '<', $fasta;
open my $out, '>', $outfile;

my %match_range;

while (my $l = <$in>) {
    chomp $l;
    my @f = split "\t", $l;
    if (@f) { # check for blank lines in input
	next if exists $match_range{$f[0]};
	$match_range{$f[0]} = join "|", $f[6], $f[7];
    }
}
close $in;

my ($scrSeqCt, $validscrSeqCt) = (0, 0);

{
    local $/ = '>';
    
    while (my $line = <$fas>) {
	$line =~ s/>//g;
	next if !length($line);
	my ($seqid, @seqs) = split /\n/, $line; 
	my $seq = join '', @seqs;
	my $useq = uc($seq);
	$scrSeqCt++ if defined $seq;
	$seqid =~ s/\s.*//;
	if (exists $match_range{$seqid}) {
	    my ($match_start, $match_end) = split /\|/, $match_range{$seqid};
	    if (defined $match_start && defined $match_end) {
		my $match_length = $match_end - $match_start;
		if ($match_length >= $length) {
		    $validscrSeqCt++;
		    my $seq_match = substr $useq, $match_start, $match_length;
		    say $out join "\n", ">".$seqid, $seq_match;
		}
	    }
	}
    }
}
close $fas;
close $out;

sub usage {
    my $script = basename($0);
  print STDERR <<END
USAGE: $script -i seqsin.fas -o filterseqs.fas -b blastfile -l length

Required:
    -i|infile       :    Fasta file of reads or contigs to filter.
    -b|blastfile    :    Blast file of query reads to some screening database.
    -o|outfile      :    A file to place the filtered sequences.

Options:
    -l|length       :    Length (integer) to be used as the lower threshold for filtering (Default: 50).
    -h|help         :    Print usage statement.
    -m|man          :    Print full documentation.
END
}
