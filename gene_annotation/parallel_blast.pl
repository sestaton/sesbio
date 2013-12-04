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

Parallel::ForkManager is a non-core Perl library that must
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

The file of sequences to BLAST. The format may be Fasta or Fastq,
and may be compressed with either gzip or bzip2.

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
'7' with is "blastxml" (BLAST XML output), and '0' which is the defaout pairwise output.

=item -e, --evalue

The e-value threshold for hits to each query. Default is 1e-5.

=item -w, --warn

Print the BLAST warnings. Defaust is no;

=item -h, --help

Print a usage statement. 

=item -m, --man

Print the full documentation.

=cut      

use 5.010;
use strict;
use warnings;
use Cwd;
use Getopt::Long;
use Pod::Usage;
use Time::HiRes qw(gettimeofday);
use File::Basename;
use File::Temp;
use IPC::System::Simple qw(system);
use Try::Tiny;
use Parallel::ForkManager;

#
# Vars with scope
#
my $infile;
my $outfile;
my $database;
my $numseqs;
my $thread;
my $cpu;
my $warn;
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
	   'w|warn'             => \$warn,
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
$blast_program //= 'blastp';
$blast_format //= 8;
$num_alignments //= 250;
$num_descriptions //= 500;
$evalue //= 1e-5;

my ($seq_files,$seqct) = split_reads($infile,$outfile,$numseqs);

open my $out, '>>', $outfile or die "\nERROR: Could not open file: $outfile\n"; 

my $pm = Parallel::ForkManager->new($thread);
$pm->run_on_finish( sub { my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data_ref) = @_;
			  for my $bl (sort keys %$data_ref) {
			      open my $report, '<', $bl or die "\nERROR: Could not open file: $bl\n";
			      print $out $_ while <$report>;
			      print $out "\n\n" if $blast_format == 0;
			      close $report;
			      unlink $bl;
			  }
			  my $t1 = gettimeofday();
			  my $elapsed = $t1 - $t0;
			  my $time = sprintf("%.2f",$elapsed/60);
			  say basename($ident)," just finished with PID $pid and exit code: $exit_code in $time minutes";
		      } );

for my $seqs (@$seq_files) {
    $pm->start($seqs) and next;
    my $blast_out = run_blast($seqs,$database,$cpu,$blast_program,
			      $blast_format,$num_alignments,
			      $num_descriptions,$evalue,$warn);
    $blasts{$blast_out} = 1;
    
    unlink $seqs;
    $pm->finish(0, \%blasts);
}

$pm->wait_all_children;

close $out;

my $t2 = gettimeofday();
my $total_elapsed = $t2 - $t0;
my $final_time = sprintf("%.2f",$total_elapsed/60);

say "\n========> Finihsed running BLAST on $seqct sequences in $final_time minutes";

exit;
#
# Subs
#
sub run_blast {
    my ($subseq_file,$database,$cpu,$blast_program,
	$blast_format,$num_alignments,$num_descriptions,
	$evalue,$warn) = @_;

    my ($dbfile,$dbdir,$dbext) = fileparse($database, qr/\.[^.]*/);
    my ($subfile,$subdir,$subext) = fileparse($subseq_file, qr/\.[^.]*/);

    my $suffix;
    if ($blast_format == 8) {
	$suffix = ".bln";
    }
    elsif ($blast_format == 7) {
	$suffix = ".blastxml";
    }
    elsif ($blast_format == 0) {
	$suffix = ".$blast_program";
    }
    my $subseq_out = $subfile."_".$dbfile.$suffix;

    my ($niceload, $blast_cmd, $exit_value);
    #$niceload  = "niceload --noswap --hard --run-mem 10g";
    #$blast_cmd = "$niceload  ".
    $blast_cmd = "blastall ". 
                 "-p $blast_program ".
	         "-e $evalue ". 
	         "-F F ". #'m S' ".      # filter simple repeats with 'seg' by default (DUST for nuc) -- Can't set w/o knowing blast program
	         "-v $num_alignments ".
	         "-b $num_descriptions ".
	         "-i $subseq_file ".
	         "-d $database ".
	         "-o $subseq_out ".
	         "-a $cpu ".
	         "-m $blast_format";
	         #"2>&1 >/dev/null'";       # we really don't want multiple processes complaining silmutaneously

    try {
	$exit_value = system([0..5],$blast_cmd);
    }
    catch {
	"\nERROR: BLAST exited with exit value $exit_value. Here is the exception: $_\n";
    };

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
    open $out, '>', $fname or die "\nERROR: Could not open file: $fname\n";
    
    push @split_files, $fname;
    my $in = get_fh($infile);
    my @aux = undef;
    my ($name, $comm, $seq, $qual);
    while (($name, $comm, $seq, $qual) = readfq(\*$in, \@aux)) {
	if ($count % $numseqs == 0 && $count > 0) {
	    $fcount++;
            $tmpiname = $iname."_".$fcount."_XXXX";
            my $fname = File::Temp->new( TEMPLATE => $tmpiname,
					 DIR => $cwd,
					 SUFFIX => ".fasta",
					 UNLINK => 0);
	    open $out, '>', $fname or die "\nERROR: Could not open file: $fname\n";

	    push @split_files, $fname;
	}
	say $out join "\n", ">".$name, $seq;
	$count++;
    }
    close $in; close $out;
    return (\@split_files, $count);
}

sub readfq {
    my ($fh, $aux) = @_;
    @$aux = [undef, 0] if (!@$aux);
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
    my ($name, $comm);
    defined $_ && do {
        ($name, $comm) = /^.(\S+)(?:\s+)(\S+)/ ? ($1, $2) : 
                	 /^.(\S+)/ ? ($1, '') : ('', '');
    };
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
    return ($name, $comm, $seq) if ($c ne '+');
    my $qual = '';
    while (<$fh>) {
        chomp;
        $qual .= $_;
        if (length($qual) >= length($seq)) {
            $aux->[0] = undef;
            return ($name, $comm, $seq, $qual);
        }
    }
    $aux->[1] = 1;
    return ($name, $seq);
}

sub get_fh {
    my ($file) = @_;

    my $fh;
    if ($file =~ /\.gz$/) {
        open $fh, '-|', 'zcat', $file or die "\nERROR: Could not open file: $file\n";
    }
    elsif ($file =~ /\.bz2$/) {
        open $fh, '-|', 'bzcat', $file or die "\nERROR: Could not open file: $file\n";
    }
    else {
        open $fh, '<', $file or die "\nERROR: Could not open file: $file\n";
    }

    return $fh;
}

sub usage {
    my $script = basename($0);
    print STDERR <<END

USAGE: $script -i seqs.fas -d db -o blast_result -n num [-t] [-a] [-b] [-v] [-p] [-bf] [-e] [-h] [-m]

Required:
    -i|infile        :    FastA/Q file to search (maybe compressed with gzip or bzip2).
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
