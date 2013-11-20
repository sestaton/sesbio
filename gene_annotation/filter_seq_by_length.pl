#!/usr/bin/env perl

=head1 NAME 
                                                                       
filter_seq_by_length.pl - Partition a set of reads for contigs by length.

=head1 SYNOPSIS    

filter_seq_by_length.pl -i seqs.fas -o filterseqs.fas -l 100 [--over]

=head1 DESCRIPTION
                                                                   
Takes as input a file of sequences and returns only those
sequences that are over (or under) a specified threshold.                                                         
                                                                       
=head1 AUTHOR 

S. Evan Staton                                                

=head1 CONTACT
 
statonse at gmail dot com

=head1 REQUIRED ARGUMENTS

=over 2

=item -i, --infile

The FastA/Q file containing reads for contigs to filter.

=item -o, --outfile

A file to place the filtered reads for contigs.

=item -l, --length

An integer specifying the length to be used as a lower
threshold for filtering.

=item --over

This option specifies that all contigs above the chosen
length will be written.

=item --under

This option specifies that all contigs under the chosen
length will be written.

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
use File::Basename;
use Pod::Usage;
       
#
# Vars 
#
my $infile;       
my $outfile;
my $length;
my $over;
my $under;
my $help;
my $man;

# counters 
my $overTotal  = 0;
my $overCount  = 0;
my $underTotal = 0;
my $underCount = 0;

GetOptions(
	   'i|infile=s'     => \$infile,
	   'o|outfile=s'    => \$outfile,
	   'l|length=i'     => \$length,
	   'over'           => \$over,
	   'under'          => \$under,
	   'h|help'         => \$help,
           'm|man'          => \$man,
    ) or pod2usage( "Try 'basename($0) --man' for more information." );

#
# Check @ARGV
#
usage() and exit(0) if $help;

pod2usage( -verbose => 2 ) if $man;

if (!$infile || !$outfile || !$length) {
    say "\nERROR: No input was given.";
    usage();      
    exit(1);
}

if ($over && $under) {
    say "\nERROR: Cannot choose both --over and --under. Exiting.";
    usage();
    exit(1);
}

# create SeqIO objects to read in and to write outfiles
my ($name, $comm, $seq, $qual);
my @aux = undef;

my $fh = get_fh($infile);
open my $out, '>', $outfile or die "\nERROR: Could not open file: $outfile";

while (($name, $comm, $seq, $qual) = readfq(\*$fh, \@aux)) {
    if ($over) {
	if ( length($seq) >= $length ) {     
	    $overCount++;
	    $overTotal += length($seq);
	    say $out join "\n", ">".$name, $seq if !defined $qual && !defined $comm;
	    say $out join "\n", ">".$name.q{ }.$comm, $seq if !defined $qual && defined $comm;
	    say $out join "\n", "@".$name, $seq, '+', $qual if defined $qual && !defined $comm;
	    say $out join "\n", "@".$name.q{ }.$comm, $seq, '+', $qual if defined $qual && defined $comm; 
	} else {	
	    $underCount++;
	    $underTotal += length($seq);
	}
    }
    if ($under) {
	if ( length($seq) < $length ) {
	    $underCount++;
            $underTotal += $seq->length;
	    say $out join "\n", ">".$name, $seq if !defined $qual && !defined $comm;
	    say $out join "\n", ">".$name.q{ }.$comm, $seq if !defined $qual && defined $comm;
	    say $out join "\n", "@".$name, $seq, '+', $qual if defined $qual && !defined $comm;
	    say $out join "\n", "@".$name.q{ }.$comm, $seq, '+', $qual if defined $qual && defined $comm;
        } else {
	    $overCount++;
            $overTotal += length($seq);
        }
    }
} 
close $fh;
close $out;

my $count = $overCount + $underCount;
my $total = $overTotal + $underTotal;
my $mean  = sprintf("%.2f", $total/$count);

say "=======================  $infile length summary";
say "Total -------> $count; Mean length $mean bp";
say "=========================================================";

if ($overCount >= 1) {
    my $overMean = sprintf("%.2f", $overTotal/$overCount); 
    say "Total number over $length bp : $overCount; Mean length: $overMean bp";
    #say "********** Sequences written to file -----------> $outfile\n";
} else {
   say "WARNING: There are no sequences over $length";
}
if ($underCount >= 1) {
    my $underMean = sprintf("%.2f", $underTotal/$underCount); 
    say "Total number under $length bp: $underCount; Mean length: $underMean bp";
} 
#say "********** Sequences written to file -----------> $outfile";

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
USAGE: $script -i seqsin.fas -o filterseqs.fas -l length [--over] [--under]

Required:
    -i|infile    :    FastA/Q file of reads or contigs to filter.
    -o|outfile   :    File to place the filtered reads or contigs.
    -l|length    :    Length (integer) to be used as the lower
                      threshold for filtering.
    --over       :    Keep only sequences over the chosen length.
                      Not to be used with the option [--under].
    --under      :    Keep only sequences under the chosen length.		      
                      Not to be used with the option [--over].  
Options:
    -h|help      :    Print usage statement.
    -m|man       :    Print full documentation.

END
}

