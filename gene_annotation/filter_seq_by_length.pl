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

The fasta file containing reads for contigs to filter.

=item -o, --outfile

A file to place the filtered reads for contigs.

=item -e, --excluded

A file to place the reads that did not pass the threshold 
for filtering.

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
use Bio::SeqIO;
use Getopt::Long;
use File::Basename;
use Pod::Usage;
       
#
# Vars 
#
my $infile;       
my $outfile;
my $excluded;      
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
	   'e|excluded=s'   => \$excluded,
	   'l|length=i'     => \$length,
	   'over'           => \$over,
	   'under'          => \$under,
	   'h|help'         => \$help,
	   'm|man'          => \$man,
) || pod2usage( "Try 'basename($0) --man' for more information." );

#
# Check @ARGV
#
usage() and exit(0) if $help;

pod2usage( -verbose => 2 ) if $man;

if (!$infile || !$outfile 
    || !$length || !$excluded) {
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
my $seqs_in = Bio::SeqIO->new('-file' => "$infile",
			      '-format' => 'fasta',
			     );

my %seqs_out = (
                'selected' => Bio::SeqIO->new('-file' => ">$outfile",
					   '-format' => 'fasta'),
                'excluded' => Bio::SeqIO->new('-file' => ">$excluded",
					    '-format' => 'fasta'),
		);

while ( my $seq = $seqs_in->next_seq() ) {
  #partition the sequences by length
    if ($over) {
	if ( $seq->length >= $length ) {     
	    $overCount++;
	    $overTotal += $seq->length;
	    $seqs_out{'selected'}->write_seq($seq);
	} else {	
	    $underCount++;
	    $underTotal += $seq->length;
	    $seqs_out{'excluded'}->write_seq($seq);
	}
    }
    if ($under) {
	if ( $seq->length < $length ) {
	    $underCount++;
            $underTotal += $seq->length;
            $seqs_out{'selected'}->write_seq($seq);
        } else {
	    $overCount++;
            $overTotal += $seq->length;
	    $seqs_out{'excluded'}->write_seq($seq);
        }
    }
} 

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
# SUBS
#
sub usage {
  my $script = basename($0);
  print STDERR <<END
USAGE: $script -i seqsin.fas -o filterseqs.fas -e undesired.fas -l length [--over] [--under]

Required:
    -i|infile    :    Fasta file of reads or contigs to filter.
    -o|outfile   :    File to place the filtered reads or contigs.
    -e|excluded  :    File to place the reads that did not pass the
                      threshold. 
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

