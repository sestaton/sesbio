#!/usr/bin/env perl

=head1 NAME 
                                                                       
parallel_blast.pl - Run multiple BLAST threads concurrently

=head1 SYNOPSIS    
 
parallel_blast.pl -i seqs.fas -o seqs_nt.bln -t 2 -n 100000 -cpu 2

=head1 DESCRIPTION
     
This script can accelerate BLAST searches by splitting an input file and 
running BLAST on multiple subsets of sequences concurrently. The size of 
the splits to make and the number of threads to create are optional. The 
input set of sequences may be in fasta or fastq format.                                                                

=head1 DEPENDENCIES

Perl 5.10 or newer. Parallel::ForkManager is a non-core Perl library that must
be installed in order for this script to work. 

Tested with:

=over

=item *
L<Parallel::ForkManger> 0.7.9 and Perl 5.8.5 (Red Hat Enterprise Linux AS release 4 (Nahant Update 9))

=item *
L<Parallel::ForkManager> 0.7.9 and Perl 5.14.2 (Red Hat Enterprise Linux Server release 5.8 (Tikanga))

=back

=head1 AUTHOR 

S. Evan Staton                                                

=head1 CONTACT
 
statonse at gmail dot com

=head1 REQUIRED ARGUMENTS

=over 2

=item -i, --infile

The file of sequences to BLAST. The format may be fasta or fastq 
but the format must be indicated (see option -sf).

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

=back

=head1 OPTIONS

=over 2

=item -t, --threads

The number of BLAST threads to spawn. Default is 1.

=item -a, --cpu

The number of processors to use for each BLAST thread. Default is 1.

=item -b, --num_aligns

The number of alignments to keep for each query. Default is 250.

=item -v, --num_desc

The number of descriptions to keep for each hit. Default is 500.

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

use v5.10;
use strict;
use warnings;
use Cwd;
use Getopt::Long;
use Pod::Usage;
use Time::HiRes qw(gettimeofday);
use File::Basename;
use File::Temp;
use IPC::Open3;
use Parallel::ForkManager;

#
# Lexical vars
#
my $infile;
my $outfile;
my $database;
my $numseqs;
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
my %blasts;    # container for reports

GetOptions(# Required
           'i|infile=s'         => \$infile,
	   'o|outfile=s'        => \$outfile,
	   'd|database=s'       => \$database,
	   'n|numseqs=i'        => \$numseqs,
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

if (!$infile || !$outfile || !$database || !$numseqs) {
    print "\nERROR: No input was given.\n";
    usage();
    exit(1);
}

#
# Set vaules
#
my $t0 = gettimeofday();
$cpu //= 1;  
$thread //= 1;

my ($seq_files,$seqct) = split_reads($infile,$outfile,$numseqs);

open(my $out, '>>', $outfile) or die "\nERROR: Could not open file: $outfile\n"; 

my $pm = Parallel::ForkManager->new($thread);
$pm->run_on_finish( sub { my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data_ref) = @_;
			  for my $bl (sort keys %$data_ref) {
			      open(my $b, '<', $bl) or die "\nERROR: Could not open file: $bl\n";
			      print $b $_ while <$report>;
			      close($report);
			      unlink($bl);
			  }
			  my $t1 = gettimeofday();
			  my $elapsed = $t1 - $t0;
			  my $time = sprintf("%.2f",$elapsed/60);
			  if ($core_dump) {
			      print "\nWARNING: It looks like $pid produced a core dump: $core_dump\n";
			  }
			  print basename($ident)," just finished with PID $pid and exit code: $exit_code in $time minutes\n";
		      } );

for my $seqs (@$seq_files) {
    $pm->start($seqs) and next;
    my $blast_out = run_blast($seqs,$database,$cpu,$blast_program,
			      $blast_format,$num_alignments,
			      $num_descriptions,$evalue);
    $blasts{$blast_out} = 1;
    
    unlink($seqs);
    $pm->finish(0, \%blasts);
}

$pm->wait_all_children;

close($out);

my $t2 = gettimeofday();
my $total_elapsed = $t2 - $t0;
my $final_time = sprintf("%.2f",$total_elapsed/60);

print "\n========> Finished running BLAST on $seqct sequences in $final_time minutes\n";

exit;
#
# Subs
#
sub run_blast {
    my ($subseq_file,$database,$cpu,$blast_program,
	$blast_format,$num_alignments,$num_descriptions,
	$evalue) = @_;

    $blast_program //= 'blastp';           
    $blast_format //= 8;
    $num_alignments //= 250;
    $num_descriptions //= 500;
    $evalue //= 1e-5;

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

    my $blast_cmd = "blastall -p $blast_program ".
	"-e $evalue ". 
	"-F T ". #'m S' ".             # filter simple repeats with 'seg' by default (DUST for nuc)
	"-v $num_alignments ".
	"-b $num_descriptions ".
	"-i $subseq_file ".
	"-d $database ".
	"-o $subseq_out ".
	"-a $cpu ".
	"-m $blast_format";

    my $bout = $subfile."_blast.out";
    my $berr = $subfile."_blast.err";
    my $pid;
    eval { $pid = open3(undef, $bout, $berr, $blast_cmd); };
    #die "open3: $@\n" if $@; 
    if ($@) { print "\nERROR: BLAST failed. If you think this is a bug, please report it. Exiting.\n"; exit(1); }
    waitpid($pid, 0);
    if ($?) {
	die "\nERROR: Child $pid exited with status of: $?. Exiting.\n";
    }
    close($bout); close($berr);
    return $subseq_out;
}

sub split_reads {
    my ($infile,$outfile,$numseqs) = @_;

    my ($iname, $ipath, $isuffix) = fileparse($infile, qr/\.[^.]*/);
    
    my $out;
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
    open($out, '>', $fname) or die "\nERROR: Could not open file: $fname\n";
    
    push @split_files, $fname;
    open(my $in, '<', $infile) or die "\nERROR: Could not open file: $infile\n";
    my @aux = undef;
    my ($name, $seq, $qual);
    while (($name, $seq, $qual) = readfq(\*$in, \@aux)) {
	if ($count % $numseqs == 0 && $count > 0) {
	    $fcount++;
            $tmpiname = $iname."_".$fcount."_XXXX";
            my $fname = File::Temp->new( TEMPLATE => $tmpiname,
                                      DIR => $cwd,
				      SUFFIX => ".fasta",
				      UNLINK => 0);
	    open($out, '>', $fname) or die "\nERROR: Could not open file: $fname\n";

	    push @split_files, $fname;
	}
	print $out join "\n", ">".$name, $seq;
	$count++;
    }
    close($in); close($out);
    return (\@split_files,$count);
}

sub readfq {
    my ($fh, $aux) = @_;
    @$aux = [undef, 0] if (!defined(@$aux));
    return if ($aux->[1]);
    if (!defined($aux->[0])) {
        while (<$fh>) {
            chomp;
            if (substr($_, 0, 1) eq '>' || substr($_, 0, 1) eq '@') {
                $aux->[0] = $_;
                last;
            }
        }
        if (!defined($aux->[0])) {
            $aux->[1] = 1;
            return;
        }
    }
    my $name;                       # This keeps paired-end information
    if (/^.?((\S+)\s(\d)\S+)/) {    # Illumina 1.8+
        $name = $1."/".$2;
    }
    elsif (/^.?(\S+)/) {            # Illumina 1.3+
        $name = $1;
    } else {
        $name = '';                 # ?
    }
    my $seq = '';
    my $c;
    $aux->[0] = undef;
    while (<$fh>) {
        chomp;
        $c = substr($_, 0, 1);
        last if ($c eq '>' || $c eq '@' || $c eq '+');
        $seq .= $_;
    }
    $aux->[0] = $_;
    $aux->[1] = 1 if (!defined($aux->[0]));
    return ($name, $seq) if ($c ne '+');
    my $qual = '';
    while (<$fh>) {
        chomp;
        $qual .= $_;
        if (length($qual) >= length($seq)) {
            $aux->[0] = undef;
            return ($name, $seq, $qual);
        }
    }
    $aux->[1] = 1;
    return ($name, $seq);
}

sub usage {
    my $script = basename($0);
    print STDERR <<END

USAGE: $script -i seqs.fas -d db -o blast_result -n num [-t] [-a] [-b] [-v] [-p] [-bf] [-e] [-h] [-m]

Required:
    -i|infile        :    Fasta file to search (contig or chromosome).
    -o|outfile       :    File name to write the blast results to.
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
