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
my $fas = get_fh($fasta);
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

my @aux = undef;
my ($name, $comm, $seq, $qual);
while (($name, $comm, $seq, $qual) = readfq(\*$fas, \@aux)) {
    $scrSeqCt++ if defined $seq;
    if (exists $match_range{$name}) {
	my ($match_start, $match_end) = split /\|/, $match_range{$name};
	if (defined $match_start && defined $match_end) {
	    my $match_length = $match_end - $match_start;
	    if ($match_length >= $length) {
		$validscrSeqCt++;
		my $seq_match = substr $seq, $match_start, $match_length;
		say $out join "\n", ">".$name, $seq_match;
	    }
	}
    }
}
close $fas;
close $out;

say STDERR "$scrSeqCt total sequences matched the target.";
say STDERR "$validscrSeqCt were above the length threshold and were written to $outfile.";

#
# methods
#
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
