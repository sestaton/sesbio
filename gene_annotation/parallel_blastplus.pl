#!/usr/bin/env perl

## NB: BLAST+ uses a lot of memory, so use this with caution. 
##     It is preferred to use the parallel_blast.pl script.

use 5.012;
use strict;
use warnings;
use Cwd;
use Getopt::Long;
use Pod::Usage;
use Time::HiRes qw(gettimeofday);
use File::Basename;
use File::Temp;
use Parallel::ForkManager;
use IPC::System::Simple qw(system);
use Try::Tiny;

#
# Vars with scope
#
my $infile;
my $outfile;
my $database;
my $numseqs;
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
    say "\nERROR: No input was given.";
    usage();
    exit(1);
}

#
# Set vaules
#
my $t0 = gettimeofday();
$cpu //= 1;
$thread //= 1;

my ($seq_files, $seqct) = split_reads($infile, $outfile, $numseqs);

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
    my $blast_out = run_blast($seqs, $database, $cpu, $blast_program, $blast_format,
			      $num_alignments, $num_descriptions, $evalue);
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
    my ($subseq_file, $database, $cpu, $blast_program, $blast_format,
	$num_alignments, $num_descriptions, $evalue) = @_;

    $blast_program //= 'blastp';
    $blast_format  //= '6';
    $num_alignments //= '250';
    $num_descriptions //= '500';
    $evalue //= '1e-5';

    my ($dbfile, $dbdir, $dbext)    = fileparse($database, qr/\.[^.]*/);
    my ($subfile, $subdir, $subext) = fileparse($subseq_file, qr/\.[^.]*/);

    my $suffix;
    if ($blast_format == 6) {
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
		    "-outfmt 6". #\"6 qseqid qlen qstart qend sseqid slen sstart send pident bitscore evalue qframe\" ".
		    "-out $subseq_out ".
		    "-num_threads $cpu ".
		    "-max_target_seqs 100000";

    my $exit_value;
    try {
        $exit_value = system([0..5],$blast_cmd);
    }
    catch {
        "\nERROR: BLAST exited with exit value $exit_value. Here is the exception: $_\n";
    };

    return $subseq_out;
}

sub split_reads {
    my ($infile, $outfile, $numseqs) = @_;

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
    -i|infile        :    FASTA/Q file to search (contig, chromosome, or reads).
    -o|outfile       :    File name to write the blast results to.
    -d|database      :    Database to search.
    -n|numseqs       :    The number of sequences to write to each split.

Options:
    -t|threads       :    Number of threads to create (Default: 1).
    -a|cpu           :    Number of processors to use for each thread (Default: 1).
    -b|num_aligns    :    Number of alignments to keep (Default: 250).
    -v|num_desc      :    Number of descriptions to keep (Default: 500).
    -p|blast_prog    :    BLAST program to execute (Default: blastp).
    -bf|blast_format :    BLAST output format (Default: 6. Type --man for more details).
    -e|evalue        :    The e-value threshold (Default: 1e-5).
    -h|help          :    Print a usage statement.
    -m|man           :    Print the full documentation.

END
}
