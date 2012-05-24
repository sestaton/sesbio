#!/usr/bin/perl -w

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

use strict;
use lib qw(/usr/local/bioperl/latest/);
use Bio::SeqIO;
use Getopt::Long;
use Pod::Usage;
use Time::HiRes qw(gettimeofday);
use File::Basename;
use File::Copy;

#
# Vars
#
my $infile;
my $num;
my $help;
my $man;

GetOptions(# Required
           'i|infile=s'   => \$infile,
	   'n|numseqs=i'  => \$num,
	  
           # Options
           'h|help'         => \$help,
           'm|man'          => \$man,
           ) || pod2usage( "Try 'basename($0) --man' for more information." );

pod2usage( -verbose => 2 ) if $man;

if (!$infile || !$num || $help) {
    print "\nERROR: No input was given.\n";
    &usage;
    exit(1);
}

# counters
my $seqct = 0;
my $seqover = 0;
my $t0 = gettimeofday();

# create SeqIO objects to read in and to write outfiles
my $seq_in  = Bio::SeqIO->new( -format => 'fasta', 
			       -file => $infile); 

my $outfile1 = $infile;
$outfile1 =~ s/\.fa.*//;
$outfile1 .= "_".$num.".fasta";
my $seqs_out = Bio::SeqIO->new( -format => 'fasta',
				-file => ">$outfile1");

my $outfile2 = $outfile1;
my $rem = "after"; 
$outfile2 =~ s/$num/$rem/;
my $seqs_out2 = Bio::SeqIO->new( -format => 'fasta',
				 -file => ">$outfile2");

while( my $seqs = $seq_in->next_seq() ) {
    if ($seqct < $num) {
	$seqct++;
	$seqs_out->write_seq($seqs);
    } else {
	$seqover++;
	$seqs_out2->write_seq($seqs);
    }
}

my $seqtot = $seqct + $seqover;
my $outfafter = $outfile1;
$outfafter =~ s/$num/$seqover/;
move("$outfile2","$outfafter");

my $t1 = gettimeofday();
my $elapsed = $t1 - $t0;
my $time = sprintf("%.2f",$elapsed/60);

print "\n$seqtot reads ($seqct in file $outfile1, $seqover in file $outfafter) written in $time minutes.\n\n";

exit;

# subs
sub usage {
    my $script = basename($0);
  print STDERR <<END
USAGE: $script -i s_1_sequence.fasta -n 100000 

Required:
    -i|infile   :    Fasta file of left paired reads.
    -n|num      :    The number of reads to select.

Options:
    -h|help     :    Print usage statement.
    -m|man      :    Print full documentation.
END
}
    
