#!/usr/bin/env perl

=head1 NAME 
                                                                       
 illumina_to_sanger.pl - Converts Illumina Phred+64 to Sanger Phred+33 

=head1 SYNOPSIS    

 illumina_to_sanger.pl -i s_1_sequence.fastq -o s_1_sequence_sanger.fastq 

=head1 DESCRIPTION
                                                                   
 Takes an Illumina fastq file as input and returns a Sanger fastq.
 This is a minimal script to just do the conversion, it gives no
 statistics or information about the files.

=head1 DEPENDENCIES

 This script uses EMBOSS and BioPerl so, both must be installed. 
 EMBOSS v6.2+ must be installed (the latest is v6.4 as of this writing), 
 as well as BioPerl v1.6.1+ (the latest is v1.6.9 as of this writing).

=head1 LICENSE

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.

=head1 AUTHOR 

S. Evan Staton                                                

=head1 CONTACT
 
statonse at gmail dot com

=head1 REQUIRED ARGUMENTS

=over 2

=item -i, --infile

The fastq file from an Illumina instrument with Phred+64 quality encoding.

=item -o, --outfile

A file to place the converted Sanger fastq reads with Phred+33 quality encoding.

=back

=head1 OPTIONS

=over 2

=item -h, --help

Print a usage statement. 

=item -m, --man

Print the full documentation.

=cut      

#
# INCLUDES    
#
use strict;
use warnings;
use Getopt::Long;
use Bio::Factory::EMBOSS;
use Pod::Usage;

#
# VARIABLE SCOPE  
#
my $infile;
my $outfile;
my $help;
my $man;

GetOptions(# Required
           '-i|infile=s'     => \$infile,
           '-o|outfile=s'    => \$outfile,
           # Options
	   '-h|help'         => \$help,
           '-m|man'          => \$man,
	   );

pod2usage( -verbose => 1 ) if $help;
pod2usage( -verbose => 2 ) if $man;

#
# CHECK @ARGV
#
if (!$infile || !$outfile) {
    print "\nERROR: No input was given.\n";
    pod2usage( -verbose => 1 );
}
    
#                       
# USE EMBOSS FOR CONVERSION
#
my $factory = Bio::Factory::EMBOSS->new;
my $seqret = $factory->program('seqret'); 
# $seqret is a Bio::Tools::Run::EMBOSSApplication object
$seqret->run({-sequence => $infile,
	      -sformat1 => 'fastq-illumina',
	      -outseq   => $outfile,
	      -osformat => 'fastq-sanger'});

exit;
