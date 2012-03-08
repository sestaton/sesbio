#!/usr/bin/perl -w

=head1 NAME 
                                                                       
pickNreads.pl - Extract N reads from a Fasta file

=head1 SYNOPSIS    

pickNreads.pl -i seqs.fasta -n 1000

=head1 DESCRIPTION
                                                                   
Takes a Fasta file and picks the first N, where N can be any number, sequences. The selection of sequences
is sequential and is based on the order in the file so, the sequences are NOT picked randomly. There are many 
programs for randomly picking sequences from a file.

The output is a single file that is named based on the selected number, N, of sequences. If the --write_all 
option is specified at the command line, a second file with n-N sequences will be written, where n is the total
number of sequences in the input.

For example if "seqs.fasta" contains 100 sequences and we want 10 we would do:

perl pickNreads.pl -i seqs.fasta -n 10

and the output would be:  seqs_10.fasta

If we want to create a split after N, and write all the sequences we would do:

perl pickNreads.pl -i seqs.fasta -n 10 --write_all

and the output would be:  seqs_10.fasta
                          seqs_90.fasta                          

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
my $write_all;
my $help;
my $man;

my $outfile1;
my $outfile2;
my $seqs_out;
my $seqs_out_over;

GetOptions(# Required
           'i|infile=s'   => \$infile,
	   'n|numseqs=i'  => \$num,
	  
           # Options
	   'write_all'    => \$write_all,
           'h|help'       => \$help,
           'm|man'        => \$man,
           ) || pod2usage( "Try 'basename($0) --man' for more information." );

pod2usage( -verbose => 2 ) if $man;

if ($help) { 
    &usage(); 
    exit(0); 
}

if (!$infile || !$num) {
    print "\nERROR: No input was given.\n";
    &usage;
    exit(1);
}

# counters
my $seq_ct = 0;
my $seq_over = 0;
my $t0 = gettimeofday();

# create SeqIO objects to read in and to write outfiles
my $seq_in  = Bio::SeqIO->new( -format => 'fasta', 
			       -file => $infile); 

$outfile1 = $infile;
$outfile1 =~ s/\.fa.*//;
$outfile1 .= "_".$num.".fasta";
$seqs_out = Bio::SeqIO->new(-format => 'fasta',
			    -file => ">$outfile1");

if ($write_all) {
    $outfile2 = $outfile1;
    my $rem = "after";
    $outfile2 =~ s/$num/$rem/;
    $seqs_out_over = Bio::SeqIO->new(-format => 'fasta',
				     -file => ">$outfile2");
}

while( my $seqs = $seq_in->next_seq() ) {
    $seq_ct++;
    if ($seq_ct <= $num) {
	$seqs_out->write_seq($seqs);
    } else {
	if ($write_all) {
	    $seq_over++;
	    $seqs_out_over->write_seq($seqs);
	}
    }
}

my $out_over = $outfile1;
$out_over =~ s/$num/$seq_over/;
move("$outfile2","$out_over") or die "Copy failed: $!" if $write_all;

my $t1 = gettimeofday();
my $elapsed = $t1 - $t0;
my $time = sprintf("%.2f",$elapsed/60);

print "\n$seq_ct reads ($num in file $outfile1, $seq_over in file $out_over) written in $time minutes.\n\n" if $write_all;
print "\n$num reads in file $outfile1 written in $time minutes.\n\n" if !$write_all;

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
    --write_all :    Write all the sequences minus [--num] to a separate file.
    -h|help     :    Print usage statement.
    -m|man      :    Print full documentation.
END
}
    
