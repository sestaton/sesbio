#!/usr/bin/perl -w

=head1 NAME 
                                                                       
parallel_blast.pl - Run multiple BLAST threads concurrently

=head1 SYNOPSIS    
 
parallel_blast.pl -i seqs.fas -o seqs_nt.bln -f fasta -t 2 -n 100000 -cpu 2

=head1 DESCRIPTION
     
Large BLAST jobs can take a very long time to complete (weeks or months
for large datasets and a large number of sequences). This script can
accelerate BLAST searches by splitting an input file and running BLAST
on multiple subsets of sequences concurrently. The size of the splits to
make and the number of threads to create are optional. The input set of
sequences may be in fasta or fastq format.                                                                

=head1 DEPENDENCIES

BioPerl and Parallel::ForkManager are non-core Perl libraries that must
be installed in order for this script to work. 

Tested with:

=over

=item *
L<BioPerl> 1.069, L<Parallel::ForkManger> 0.7.9 and Perl 5.8.5 (Red Hat Enterprise Linux AS release 4 (Nahant Update 9))

=back

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

The file of sequences to BLAST. The format may be fasta or fastq 
but the format must be indicated (see option -f).

=item -o, --outfile

A file to place the BLAST results.

=item -n, --numseqs

The size of the splits to create. This number determines how many 
sequences will be written to each split. 

NB: If the input sequence file has millions of sequences and a 
very small number is given fo the split value then there could 
potentially be hundreds of thousands of files created. 

=item -d, --database

The BLAST database to search. 

=item -sf, --seq_format

The format of the input sequences. Must be one of fasta or fastq.

=back

=head1 OPTIONS

=over 2

=item -t, --threads

The number of BLAST threads to spawn. Default is 1.

=item -a, --cpu

The number of processors to use for each BLAST thread. Default is 1.

=item -b, --num_aligns

The number of alignments to keep for each query. Default is 100000.

=item -v, --num_desc

The number of descriptions to keep for each hit. Default is 100000.

=item -p, --blast_prog

The BLAST program to execute. Default is blastp.

=item -bf, --blast_format

The BLAST output format. Default is 8.
NB: The only allowed options are '8' which is "blasttable" (tabular BLAST output),
or '7' with is "blastxml" (BLAST XML output).

=item -e, --evalue

The e-value threshold for hits to each query. Default is 1e-5.

=item -h, --help

Print a usage statement. 

=item -m, --man

Print the full documentation.

=cut      

# TODO: Add BLAST options
#       Add BLAST progress meter       


use strict;
use Getopt::Long;
use Pod::Usage;
use Time::HiRes qw(gettimeofday);
use File::Basename;
use File::Temp;
use File::Spec;
use lib qw(/iob_home/jmblab/statonse/apps/perlmod/Parallel-ForkManager-0.7.9/lib);
use Parallel::ForkManager;
use Bio::SeqIO;
use POSIX qw(ceil);
use List::MoreUtils qw(uniq);

#
# Vars with scope
#
my $infile;
my $outfile;
my $database;
my $numseqs;
my $format;
my $thread;
my $cpu;
my $help;
my $man;

# BLAST-specific variables
my $blast_program;
my $blast_format;
my $num_alignments;
my $num_descriptions;
my $evalue;

GetOptions(# Required
           "i|infile=s"         => \$infile,
	   "o|outfile=s"        => \$outfile,
	   "d|database=s"       => \$database,
	   "n|numseqs=i"        => \$numseqs,
	   "sf|seq_format=s"    => \$format,
	   # Options
	   "a|cpu=i"            => \$cpu,
	   "b|num_aligns=i"     => \$num_alignments,
	   "v|num_desc=i"       => \$num_descriptions,
	   "p|blast_prog=s"     => \$blast_program,
	   "e|evalue=f"         => \$evalue,
	   "bf|blast_format=i"  => \$blast_format,
	   "t|threads=i"        => \$thread,
           "h|help"             => \$help,
           "m|man"              => \$man,
           ) || pod2usage( "Try '$0 --man' for more information." );

if ($help) {
    &usage();
    exit(0);
};

pod2usage( -verbose => 2 ) if $man;

if (!$infile || !$format || 
    !$outfile || !$database || 
    !$numseqs) {
    print "\nERROR: No input was given.\n";
    &usage();
    exit(0);
}

#
# Set vaules
#
my $t0 = gettimeofday();
$cpu = defined($cpu) ? $cpu : '1';          # we are going to set defaults this way
$thread = defined($thread) ? $thread : '1'; # to work with Perl versions released prior to 5.10

my ($seq_files,$seqct) = split_reads($infile,$outfile,$numseqs,$format);

my $pm = Parallel::ForkManager->new($thread);
$pm->run_on_finish( sub { my ($pid, $exit_code, $ident) = @_;
			  my $t1 = gettimeofday();
			  my $elapsed = $t1 - $t0;
			  my $time = sprintf("%.2f",$elapsed/60);
			  print "$ident just finished with PID $pid and exit code: $exit_code in $time minutes\n";
		        } );

open(my $out, '>>', $outfile) or die "\nERROR: Could not open file: $outfile\n";

foreach my $seqs (@$seq_files) {
    $pm->start($seqs) and next;
    my $blast_out = run_blast($seqs,$database,$cpu,$blast_program,$blast_format,$num_alignments,$num_descriptions,$evalue);

    open(my $report, '<', $blast_out) or die "\nERROR: Could not open file: $blast_out\n";
    while(my $line = <$report>) {
	print $out $line;
    }
    close($report);
    unlink($blast_out);

    unlink($seqs);
    $pm->finish;
}
#printf "%s: Waiting for some child to finish\n", scalar localtime;
$pm->wait_all_children;

#printf "%s: All processes finished.\n", scalar localtime;

close($out);

my $t2 = gettimeofday();
my $total_elapsed = $t2 - $t0;
my $final_time = sprintf("%.2f",$total_elapsed/60);

print "\n========> Finihsed running BLAST on $seqct sequences in $final_time minutes\n";

exit;
#
# Subs
#
sub run_blast {
    
    my ($subseq_file,$database,$cpu,$blast_program,$blast_format,$num_alignments,$num_descriptions,$evalue) = @_;

    $blast_program = defined($blast_program) ? $blast_program : 'blastp';          # we can set defaults with much less typing 
    $blast_format  = defined($blast_format) ? $blast_format : '8';                 # if we 'use 5.010' but we'll try to be
    $num_alignments = defined($num_alignments) ? $num_alignments : '100000';       # compatible
    $num_descriptions = defined($num_descriptions) ? $num_descriptions : '100000';
    $evalue = defined($evalue) ? $evalue : '1e-5';

    my ($dbfile,$dbdir,$dbext) = fileparse($database, qr/\.[^.]*/);

    my $subseq_out = $subseq_file;
    $subseq_out .= "_".$dbfile.".bln";

    my $blast_cmd = "blastall -p $blast_program ".
	            "-e $evalue ". 
		    "-F 'm S' ".             # filter simple repeats with 'seg' by default (DUST for nuc)
		    "-v $num_alignments ".
		    "-b $num_descriptions ".
		    "-i $subseq_file ".
		    "-d $database ".
		    "-o $subseq_out ".
		    "-a $cpu ".
		    "-m $blast_format";

    system($blast_cmd);
    return($subseq_out);

}

sub split_reads {

    my ($input,$output,$numseqs,$format) = @_;

    my ($iname, $ipath, $isuffix) = fileparse($input, qr/\.[^.]*/);
    
    my $seq_in  = Bio::SeqIO->new(-file  => $input,
				  -format => $format);
    my $count = 0;
    my $fcount = 1;
    my @split_files;
    $iname =~ s/\.fa.*//;     # clean up file name like seqs.fasta.1
    
    my $cwd = File::Spec->curdir();

    my $tmpiname = $iname."_".$fcount."_XXXX";
    my $fname = File::Temp->new( TEMPLATE => $tmpiname,
                                 DIR => $cwd,
				 UNLINK => 0);
    
    my $seq_out = Bio::SeqIO->new(-file => ">$fname", 
      			          -format=>'fasta');

    push(@split_files,$fname);
    while (my $seq = $seq_in->next_seq) {
	if ($count % $numseqs == 0 && $count > 0) {
	    $fcount++;
            $tmpiname = $iname."_".$fcount."_XXXX";
            $fname = File::Temp->new( TEMPLATE => $tmpiname,
                                      DIR => $cwd,
				      UNLINK => 0);

	    $seq_out = Bio::SeqIO->new(-file => ">$fname", 
				       -format=>'fasta');

	    push(@split_files,$fname);
	}
	$seq_out->write_seq($seq);
	$count++;
    }

    my $filect = ceil($count/$numseqs);

    my @unique_split_files = uniq(@split_files);
    return (\@unique_split_files,$count);
}

sub usage {
    my $script = basename($0);
    print STDERR <<END

USAGE: $script -i seqs.fas -d db -sf fasta|fastq -o blast_result -n num [-t] [-a] [-b] [-v] [-p] [-bf] [-e]

Required:
    -i|infile        :    Fasta file to search (contig or chromosome).
    -o|outfile       :    File name to write the blast results to.
    -sf|seq_format   :    The format of the sequence (must be "fasta" or "fastq").
    -d|database      :    Database to search.
    -n|numseqs       :    The number of sequences to write to each split.

Options:
    -t|threads       :    Number of threads to create (Default: 1).
    -a|cpu           :    Number of processors to use for each thread (Default: 1).
    -b|num_aligns    :    Number of alignments to keep (Default: 100000).
    -v|num_desc      :    Number of descriptions to keep (Default: 100000).
    -p|blast_prog    :    BLAST program to execute (Default: blastp).
    -bf|blast_format :    BLAST output format (Default: 8. Type --man for more details).
    -e|evalue        :    The e-value threshold (Default: 1e-5).

END
}
