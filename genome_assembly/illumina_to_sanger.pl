#!/usr/bin/env perl

=head1 NAME 
                                                                       
 illumina_to_sanger.pl - Converts Illumina Phred+64 to Sanger Phred+33 

=head1 SYNOPSIS    

 illumina_to_sanger.pl -i s_1_sequence.fastq -o s_1_sequence_sanger.fastq 

=head1 DESCRIPTION
                                                                   
 Takes an Illumina FASTQ file as input and returns a Sanger FASTQ.
 This is a minimal script to just do the conversion, it gives no
 statistics or information about the files.

=head1 DEPENDENCIES

 This script uses EMBOSS and BioPerl so, both must be installed. 
 EMBOSS v6.2+ must be installed (the latest is v6.4 as of this writing), 
 as well as BioPerl v1.6.1+ (the latest is v1.6.9 as of this writing).

=head1 AUTHOR 

S. Evan Staton                                                

=head1 CONTACT
 
statonse at gmail dot com

=head1 REQUIRED ARGUMENTS

=over 2

=item -i, --infile

The FASTQ file from an Illumina instrument with Phred+64 quality encoding.

=item -o, --outfile

A file to place the converted Sanger FASTQ reads with Phred+33 quality encoding.

=back

=head1 OPTIONS

=over 2

=item -h, --help

Print a usage statement. 

=item -m, --man

Print the full documentation.

=cut      

use strict;
use warnings;
use Getopt::Long;
use Bio::Factory::EMBOSS;
use Pod::Usage;

my $infile;
my $outfile;
my $help;
my $man;

GetOptions(# Required
           'i|infile=s'  => \$infile,
           'o|outfile=s' => \$outfile,
           # Options
	   'h|help'      => \$help,
           'm|man'       => \$man,
	   );

pod2usage( -verbose => 1 ) if $help;
pod2usage( -verbose => 2 ) if $man;

if (!$infile || !$outfile) {
    print "\nERROR: No input was given.\n";
    pod2usage( -verbose => 1 );
}

my $factory = Bio::Factory::EMBOSS->new;
my $seqret = $factory->program('seqret'); 

# $seqret is a Bio::Tools::Run::EMBOSSApplication object
$seqret->run({-sequence => $infile,
	      -sformat1 => 'fastq-illumina',
	      -outseq   => $outfile,
	      -osformat => 'fastq-sanger'});
