#!/usr/bin/env perl

=head1 NAME 
                                                                       
pickNreads.pl - Extract N reads from a Fasta file

=head1 SYNOPSIS    

pickNreads.pl -i seqs.fasta -n 1000

=head1 DESCRIPTION
                                                                   
Takes a Fasta file and picks the first N, where N
can be any number, sequences. The selection of sequences
is sequential and is based on the order in the file so, 
the sequences are NOT picked randomly. There are many 
programs for randomly picking sequences from a file.

The output is two files where the first is named
based on the selected number of sequences and the second
will contain every other sequence in the file. For
example if "seqs.fasta" contains 100 sequences and we 
want 10 we would do:

perl pickNreads.pl -i seqs.fasa -n 10

and the output would be:  seqs_10.fasta
                          seqs_90.fasfa

=head1 AUTHOR 

S. Evan Staton                                                

=head1 CONTACT
 
statonse at gmail dot com

=head1 REQUIRED ARGUMENTS

=over 2

=item -i, --infile

The Fasta file from which to select sequences.

=item -n, --numreads

The exact number of sequences to select.

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
use Getopt::Long;
use Pod::Usage;
use File::Basename;

#
# Vars
#
my $infile;
my $outfile;
my $num;
my $help;
my $man;

GetOptions(# Required
           'i|infile=s'   => \$infile,
	   'o|outfile=s'  => \$outfile,
	   'n|numseqs=i'  => \$num,
	  
           # Options
           'h|help'         => \$help,
           'm|man'          => \$man,
           ) || pod2usage( "Try 'basename($0) --man' for more information." );

pod2usage( -verbose => 2 ) if $man;
usage() and exit(0) if $help;

if (!$infile || !$num || !$outfile) {
    print "\nERROR: No input was given.\n";
    usage();
    exit(1);
}

# counters
my $seqct = 0;

my ($name, $comm, $seq, $qual);
my @aux = undef;

my $fh = get_fh($infile);
open my $out, '<', $outfile or die "\nERROR: Could not open file: $outfile";

while (($name, $comm, $seq, $qual) = readfq(\*$fh, \@aux)) {
    if ($seqct < $num) {
	$seqct++;
	if (defined $qual) {
	    say join "\n", "@".$name, $seq, '+', $qual;
	}
	else {
	    say join "\n", ">".$name, $seq;
	}
    }
}

exit;
#
# methods
#
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
USAGE: $script -i s_1_sequence.fasta -n 100000 

Required:
    -i|infile   :    Fasta file of r
    -n|num      :    The number of reads to select.

Options:
    -h|help     :    Print usage statement.
    -m|man      :    Print full documentation.
END
}
    
