#!/usr/bin/perl -w

=head1 NAME 
                                                                       
fq2fas.pl - Convert a Fastq file to Fasta

=head1 SYNOPSIS    

fq2fas.pl -i seqs.fastq -o seqs.fasta

=head1 DESCRIPTION
                                                                   
Takes a Fastq file (any variant) and returns a Fasta file.

=head1 AUTHOR 

S. Evan Staton                                                

=head1 CONTACT
 
statonse at gmail dot com

=head1 REQUIRED ARGUMENTS

=over 2

=item -i, --infile

The Fastq file from an Illumina instrument (encoding is irrelevant).

=item -o, --outfile

A file to place the converted Fasta reads.

=back

=head1 OPTIONS

=over 2

=item -h, --help

Print a usage statement. 

=item -m, --man

Print the full documentation.

=cut   

use strict;
use Getopt::Long;
use Pod::Usage;
use Time::HiRes qw(gettimeofday);
use File::Basename;

#
# VARIABLE SCOPE  
#
my $infile;
my $outfile;
my $help;
my $man;

GetOptions(# Required
           "-i|infile=s"     => \$infile,
           "-o|outfile=s"    => \$outfile,
           # Options
           "-h|help"         => \$help,
           "-m|man"          => \$man,
           ) || pod2usage( "Try 'basename($0) --man' for more information." );

&usage if $help;
pod2usage( -verbose => 2 ) if $man;

if (!$infile || !$outfile) {
    print "\nERROR: No input was given.\n";
    &usage;
    exit(0);
}

# counters
my $seqct = 0;
my $t0 = gettimeofday();

open(my $FQ,'<',$infile) or die "\nERROR: Could not open file: $infile\n"; 
open(my $FA,'>',$outfile) or die "\nERROR: Could not open file: $outfile\n";

while(my $h = <$FQ>) {
    $seqct++;
    $h =~ s/\@/\>/;
    my $s = <$FQ>;
    my $h2 = <$FQ>;
    my $q = <$FQ>;
    print $FA $h.$s;
}

close($FQ);
close($FA);

my $t1 = gettimeofday();
my $elapsed = $t1 - $t0;
my $time = sprintf("%.2f",$elapsed/60);

print "\n$seqct reads converted to Fasta in file $outfile in $time minutes.\n\n";

exit;

#
# subs
#
sub usage {
    my $script = basename($0);
  print STDERR <<END
USAGE: $script -i s_1_sequence.fastq -o s_1_sequence.fasta 

Required:
    -i|infile       :    Fastq file of reads.
    -o|outfile      :    File to place the converted Fasta reads.

Options:
    -h|help         :    Print usage statement.
    -m|man          :    Print full documentation.
END
}
    
