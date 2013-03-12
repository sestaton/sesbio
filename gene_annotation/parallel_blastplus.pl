#!/usr/bin/env perl

use v5.12;
use strict;
use warnings;
use Cwd;
use Getopt::Long;
use Pod::Usage;
use Time::HiRes qw(gettimeofday);
use File::Basename;
use File::Temp;
use Parallel::ForkManager;
use Bio::SeqIO;

#
# Vars with scope
#
my $infile;
my $outfile;
my $database;
my $numseqs;
my $format;
my $thread;
my $max_target_seqs;
my $cpu;
my $help;
my $man;

# BLAST-specific variables
my $blast_program;
my $blast_format;
my $num_alignments;
my $num_descriptions;
my $evalue;
my %blasts;    # container for reports

GetOptions(# Required
           'i|infile=s'         => \$infile,
	   'o|outfile=s'        => \$outfile,
	   'd|database=s'       => \$database,
	   'n|numseqs=i'        => \$numseqs,
	   'sf|seq_format=s'    => \$format,
	   # Options
	   'a|cpu=i'            => \$cpu,
	   'b|num_aligns=i'     => \$num_alignments,
	   'v|num_desc=i'       => \$num_descriptions,
	   'p|blast_prog=s'     => \$blast_program,
	   'e|evalue=f'         => \$evalue,
	   'bf|blast_format=i'  => \$blast_format,
	   't|threads=i'        => \$thread,
           'h|help'             => \$help,
           'm|man'              => \$man,
           ) || pod2usage( "Try '$0 --man' for more information." );

#
# Check @ARGV
#
usage() and exit(0) if $help;

pod2usage( -verbose => 2 ) if $man;

if (!$infile || !$format || 
    !$outfile || !$database || 
    !$numseqs) {
    print "\nERROR: No input was given.\n";
    usage();
    exit(1);
}

#
# Set vaules
#
my $t0 = gettimeofday();
$cpu = defined($cpu) ? $cpu : '1';          # we are going to set defaults this way
$thread = defined($thread) ? $thread : '1'; # to work with Perl versions released prior to 5.10

my ($seq_files,$seqct) = split_reads($infile,$outfile,$numseqs,$format);

open(my $out, '>>', $outfile) or die "\nERROR: Could not open file: $outfile\n"; 

my $pm = Parallel::ForkManager->new($thread);
$pm->run_on_finish( sub { my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data_ref) = @_;
			  foreach my $bl (sort keys %$data_ref) {
			      open(my $report, '<', $bl) or die "\nERROR: Could not open file: $bl\n";
			      while(my $line = <$report>) {
				  print $out $line;
			      }
			      close($report);
			      unlink($bl);
			  }
			  my $t1 = gettimeofday();
			  my $elapsed = $t1 - $t0;
			  my $time = sprintf("%.2f",$elapsed/60);
			  print "$ident just finished with PID $pid and exit code: $exit_code in $time minutes\n";
		      } );

foreach my $seqs (@$seq_files) {
    $pm->start($seqs) and next;
    my $blast_out = run_blast($seqs,$database,$cpu,$blast_program,$blast_format,$num_alignments,$num_descriptions,$evalue);
    $blasts{$blast_out} = 1;
    
    unlink($seqs);
    $pm->finish(0, \%blasts);
}

$pm->wait_all_children;

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

    $blast_program = defined($blast_program) ? $blast_program : 'blastp';          # We can set defaults with much less typing 
    $blast_format  = defined($blast_format) ? $blast_format : '8';                 # if we 'use 5.010' but we'll try to be compatible.
    $num_alignments = defined($num_alignments) ? $num_alignments : '250';          # These are the BLAST defaults, increase as needed       
    $num_descriptions = defined($num_descriptions) ? $num_descriptions : '500';    # e.g., for OrthoMCL.
    $evalue = defined($evalue) ? $evalue : '1e-5';

    my ($dbfile,$dbdir,$dbext) = fileparse($database, qr/\.[^.]*/);
    my ($subfile,$subdir,$subext) = fileparse($subseq_file, qr/\.[^.]*/);

    my $suffix;
    if ($blast_format == 8) {
	$suffix = ".bln";
    }
    elsif ($blast_format == 7) {
	$suffix = ".blastxml";
    }
    my $subseq_out = $subfile."_".$dbfile.$suffix;

    my $blast_cmd = "blastn ".
	            "-dust no ".
	            "-query $subseq_file ".
		    "-db $database ".
		    "-outfmt \"6 qseqid qlen qstart qend sseqid slen sstart send pident bitscore evalue qframe\" ". 
		    "-out $subseq_out ".
		    "-num_threads $cpu ". 
		    "-max_target_seqs 100000";

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
    
    my $cwd = getcwd();

    my $tmpiname = $iname."_".$fcount."_XXXX";
    my $fname = File::Temp->new( TEMPLATE => $tmpiname,
                                 DIR => $cwd,
				 SUFFIX => ".fasta",
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
				      SUFFIX => ".fasta",
				      UNLINK => 0);

	        $seq_out = Bio::SeqIO->new(-file => ">$fname", 
					   -format=>'fasta');

	    push(@split_files,$fname);
	}
	$seq_out->write_seq($seq);
	$count++;
    }

    return (\@split_files,$count);
}

sub usage {
    my $script = basename($0);
    print STDERR <<END

USAGE: $script -i seqs.fas -d db -sf fasta|fastq -o blast_result -n num [-t] [-a] [-b] [-v] [-p] [-bf] [-e] [-h] [-m]

Required:
    -i|infile        :    Fasta file to search (contig or chromosome).
    -o|outfile       :    File name to write the blast results to.
    -sf|seq_format   :    The format of the sequence (must be "fasta" or "fastq").
    -d|database      :    Database to search.
    -n|numseqs       :    The number of sequences to write to each split.

Options:
    -t|threads       :    Number of threads to create (Default: 1).
    -a|cpu           :    Number of processors to use for each thread (Default: 1).
    -b|num_aligns    :    Number of alignments to keep (Default: 250).
    -v|num_desc      :    Number of descriptions to keep (Default: 500).
    -p|blast_prog    :    BLAST program to execute (Default: blastp).
    -bf|blast_format :    BLAST output format (Default: 8. Type --man for more details).
    -e|evalue        :    The e-value threshold (Default: 1e-5).
    -h|help          :    Print a usage statement.
    -m|man           :    Print the full documentation. 

END
}
