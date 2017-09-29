#!/usr/bin/env perl

=head1 NAME 
                                                                       
qseq2fq.pl - Convert a qseq file to Fastq

=head1 SYNOPSIS    

qseq2fq.pl -i s_1_1_qseq.txt -o s_1_1.fastq 

=head1 DESCRIPTION
                                                                   
Convert a qseq file from an Illumina instrument to Fastq format. Optionally,
filter out the reads that did not pass the Chastity filter.

=head1 AUTHOR 

S. Evan Staton                                                

=head1 CONTACT
 
evan at evanstaton dot com

=head1 REQUIRED ARGUMENTS

=over 2

=item -i, --infile

The qseq file to be converted.

=item -o, --outfile

A file to place the converted Illumina reads.

=back

=head1 OPTIONS

=over 2

=item -c, --chastity

Filter out the reads that did not pass the Chastity filter when
doing the conversion.

=item -d, --description

Print the full description line in the Fastq file. 

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
use Time::HiRes qw(gettimeofday);

my $qseq;
my $fastq;
my $chastity;
my $description;
my $man;
my $help;

GetOptions(
	   'i|infile=s'    => \$qseq,
	   'o|outfile=s'   => \$fastq,
           'c|chastity'    => \$chastity,
	   'd|description' => \$description,
	   'm|man'         => \$man,
	   'h|help'        => \$help,
	   );

pod2usage( -verbose => 2 ) if $man;

#
# Check @ARGV
#
usage() and exit(0) if $help;

if (!$qseq || !$fastq) {
    say "\nERROR: Command line not parsed correctly.";
    usage();
    exit(1);
}

my $t0 = gettimeofday();
my $readnum = 0;
my $filterednum = 0;

open my $qs, '<', $qseq or die "\nERROR: Could not open file: $qseq\n";
open my $fq, '>', $fastq or die "\nERROR: Could not open file: $fastq\n";

while (my $line = <$qs>) {
    chomp $line;
    my @reads = split /\t/, $line;
    $readnum++;
    if ($chastity && $description) {
	if ($reads[10] == 1) {
	    $filterednum++;
	    say $fq "@","$reads[0]:$reads[2]:$reads[3]:$reads[4]:$reads[5]#$reads[6]/$reads[7]";
	    say $fq "$reads[8]";
	    say $fq "+","$reads[0]:$reads[2]:$reads[3]:$reads[4]:$reads[5]#$reads[6]/$reads[7]";
	    say $fq "$reads[9]";
	}
    }
    if ($chastity && !$description) {
	if ($reads[10] == 1) {
	    $filterednum++;
            say $fq "@","$reads[0]:$reads[2]:$reads[3]:$reads[4]:$reads[5]#$reads[6]/$reads[7]";
            say $fq "$reads[8]";
            say $fq "+";
            say $fq "$reads[9]";
        }
    }
    if (!$chastity && $description) {
	say $fq "@","$reads[0]:$reads[2]:$reads[3]:$reads[4]:$reads[5]#$reads[6]/$reads[7]";
        say $fq "$reads[8]";
	say $fq "+","$reads[0]:$reads[2]:$reads[3]:$reads[4]:$reads[5]#$reads[6]/$reads[7]";
	say $fq "$reads[9]";
    }
    if (!$chastity && !$description) {
	say $fq "@","$reads[0]:$reads[2]:$reads[3]:$reads[4]:$reads[5]#$reads[6]/$reads[7]";
	say $fq "$reads[8]";
	say $fq "+";
	say $fq "$reads[9]";
    }
}
close $qs;
close $fq;

my $t1 = gettimeofday();
my $elapsed = $t1 - $t0;
my $time = sprintf("%.2f",$elapsed/60);

if ($chastity) {
    say "\n========== Done. $readnum reads total. $filterednum Chastity filtered reads converted in $time minutes.";
} 
else {
    say "\n========== Done. $readnum reads converted in $time minutes.";
}

exit; 
#
# Methods
#
sub usage {
    my $script = basename($0);
  print STDERR <<END
USAGE: $script -i s_1_1_qseq.txt -o s_1_1.fastq 

Required:
    -i|infile      :    Qseq file of reads to filter and convert.
    -o|outfile     :    File to place the filtered reads (in Fastq format).
    
Options:
    -c|chastity    :    Only print out the reads that passed the Chastity filter.
    -d|description :    Print out the full description line after the "+" (though some
                        like to keep this line it is really unnecessary).
    -h|help        :    Print usage statement.
    -m|man         :    Print full documentation.
END
}
