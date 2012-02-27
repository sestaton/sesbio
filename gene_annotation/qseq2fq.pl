#!/usr/bin/perl -w

=head1 NAME 
                                                                       
qseq2fq.pl - Convert a qseq file to Fastq

=head1 SYNOPSIS    

qseq2fq.pl -i s_1_1_qseq.txt -o s_1_1.fastq 

=head1 DESCRIPTION
                                                                   
Convert a qseq file from an Illumina instrument to Fastq format. Optionally,
filter out the reads that did not pass the Chastity filter.

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

#
# Includes
#
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use Time::HiRes qw(gettimeofday);

#
# Vars
#
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
# Check @ARGVs
#
if (!$qseq || !$fastq || $help) {
    print "\nERROR: Command line not parsed correctly.\n";
    &usage();
    exit(1);
}

#
# Counters
#
my $t0 = gettimeofday();
my $readnum = 0;
my $filterednum = 0;

open(my $qs, '<', $qseq) or die "\nERROR: Could not open file: $qseq\n";
open(my $fq, '>', $fastq) or die "\nERROR: Could not open file: $fastq\n";

while (my $line = <$qs>) {
    chomp $line;
    my @reads = split(/\t/,$line);
    $readnum++;
    if ($chastity && $description) {
	if ($reads[10] == 1) {
	    $filterednum++;
	    print $fq "@","$reads[0]:$reads[2]:$reads[3]:$reads[4]:$reads[5]#$reads[6]/$reads[7]\n";
	    print $fq "$reads[8]\n";
	    print $fq "+","$reads[0]:$reads[2]:$reads[3]:$reads[4]:$reads[5]#$reads[6]/$reads[7]\n";
	    print $fq "$reads[9]\n";
	}
    }
    if ($chastity && !$description) {
	if ($reads[10] == 1) {
	    $filterednum++;
            print $fq "@","$reads[0]:$reads[2]:$reads[3]:$reads[4]:$reads[5]#$reads[6]/$reads[7]\n";
            print $fq "$reads[8]\n";
            print $fq "+\n";
            print $fq "$reads[9]\n";
        }
    }
    if (!$chastity && $description) {
	print $fq "@","$reads[0]:$reads[2]:$reads[3]:$reads[4]:$reads[5]#$reads[6]/$reads[7]\n";
	print $fq "$reads[8]\n";
	print $fq "+","$reads[0]:$reads[2]:$reads[3]:$reads[4]:$reads[5]#$reads[6]/$reads[7]\n";
	print $fq "$reads[9]\n";
    }
    if (!$chastity && !$description) {
	print $fq "@","$reads[0]:$reads[2]:$reads[3]:$reads[4]:$reads[5]#$reads[6]/$reads[7]\n";
	print $fq "$reads[8]\n";
	print $fq "+\n";
	print $fq "$reads[9]\n";
    }
}
close($qs);
close($fq);

my $t1 = gettimeofday();
my $elapsed = $t1 - $t0;
my $time = sprintf("%.2f",$elapsed/60);

if ($chastity) {
    print "\n========== Done. $readnum reads total. $filterednum Chastity filtered reads converted in $time minutes.\n";
} else {
    print "\n========== Done. $readnum reads converted in $time minutes.\n";
}

exit; 

#
# Subs
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
