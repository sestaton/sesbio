#!/usr/bin/env perl

#NB: This is a rewrite of: https://github.com/sestaton/transposome-scripts/blob/master/full_analysis_with_varying_coverage.pl
# for comparing the performance of Transposome at different levels of coverage.

## pragmas and library imports
use 5.014;
use strict;
use warnings;
use autodie;
use Getopt::Long;
use File::Basename;
use File::Find;
use File::Path          qw(make_path remove_tree);
use IPC::System::Simple qw(system);
use Cwd                 qw(getcwd);
use Regexp::Common;
use Sort::Naturally;
#use Data::Dump::Color;

my %opts;

GetOptions(\%opts,
    'accession|id=s',
    'sequence_format|sf=s',  
    'output_dir|p=s',        
    'repeat_database|db=s',  
    'threads|t=i',           
    'percent_identity|pid=i',
    'coverage|cov=f',        
    'cluster_size|cls=i',    
    'species|s=s',
   'h|help',                
);

## check input
usage() and exit(0) if $opts{help};

if (!$opts{accession} || !$opts{output_dir} || !$opts{repeat_database} || !$opts{species}) {
    say "\n[ERROR]: Command line not parsed correctly. Check input and refer to the usage. Exiting.\n";
    usage();
    exit(1);
}

if (-e $opts{output_dir}) {
    say STDERR "\nERROR: $opts{output_dir} exists. Remove it or rename it and try again. Exiting.\n";
    exit(1);
}
else {
    make_path($opts{output_dir}, {verbose => 0, mode => 0711,});
}

my $outfile = File::Spec->catfile($opts{output_dir}, $opts{accession}.'_all_transposome_results.tsv');
open my $out, '>>', $outfile;
say $out join "\t", 'X-coverage (number_clustered)', 'number_unclustered', 'repeat_fraction', 
    'singleton_fraction', 'total_repeat_percentage', 'total_annotated_repeat_percentage', 
    'family_count', 'elapsed_time', 'memory_usage';

$opts{threads} //= 1;
$opts{cluster_size} //= 100;
$opts{percent_identity}  //= 90;
$opts{coverage} //= 0.55;
$opts{sequence_format} //= 'fasta';
my $merge_thresh = 100;

# 2C Mbp C-values
my %genome_sizes = ( dmel => 165, zmays => 2665, mmus => 2640, athal => 156 );
unless (exists $genome_sizes{ $opts{species} }) {
    say STDERR "\nERROR: '$opts{species} not recognized. Must be one of: dmel, zmays, mmus. Exiting.\n";
}

my ($forward, $reverse) = ($opts{accession}.'_1_p.fastq.gz', $opts{accession}.'_2_p.fastq.gz');
unless (-e $forward && -e $reverse) {
    say STDERR "\nERROR: $forward or $reverse do not exist. Exiting.\n";
    exit(1);
}

my ($iname, $ipath, $isuffix) = fileparse($forward, qr/\.[^.]*/);
my $base = $iname =~ s/_1_p.fastq//r;

# instead of simply generating N-reads (commented out), we generate
# X-coverage for each species
my %cvalues;
#for (my $i = 1e5; $i <= 2e6; $i += 2e5) {
for my $sp (keys %genome_sizes) { 
    my $cvalue = $genome_sizes{$sp} * 1e6;
    for my $x (qw(0.001 0.003 0.005 0.007 0.009 0.01 0.03 0.05)) {
	my $reads = $x*$cvalue / 100;
	my $bases = $x*$cvalue;
	my $pairs = $reads / 2;
	$cvalues{$sp}{$x} =  { bases => $bases, reads => $reads, pairs => $pairs };
    }
}

#}
#dd \%cvalues and exit;

for my $xcov (sort { $a <=> $b } keys %{$cvalues{ $opts{species} }}) {
    say STDERR "=====> Working on $xcov for $opts{species}...";
    $merge_thresh += 100;

    my $intseq = join_pairs($cvalues{ $opts{species} }{$xcov}{reads}, $forward, $reverse);

    my $outdir   = File::Spec->catdir($opts{output_dir}, $base."_${xcov}X");
    my $logfile  = $base.'_transposome_'."${xcov}X".'.log';       
    my $clogfile = $base.'_transposome_clustering_'."${xcov}X".'.log';
    # for getting family counts
    my $sumfile  = File::Spec->catfile($outdir, $base.'_transposome_clustering_'."${xcov}X".'_annotations_summary.tsv'); 

    my %run_opts = (
	coverage        => $xcov,
	sequence_file   => $intseq,
	sequence_format => $opts{sequence_format}, 
	threads         => $opts{threads},
	outdir          => $outdir,
	repeatdb        => $opts{repeat_database},
	cluster_size    => $opts{cluster_size},
	logfile         => $logfile,
	cluster_logfile => $clogfile);

    my $config = write_config(\%run_opts);
    system('valgrind', '--tool=massif', 'transposome', '-c', $config) == 0 or die $!;
    my $mem = get_mem_usage();

    my $reslog = File::Spec->catfile($outdir, $logfile);
    write_results($xcov, $cvalues{ $opts{species} }{$xcov}{reads}, $reslog, $sumfile, $mem, $out);
    remove_tree($outdir, { safe=> 1 });
    unlink $intseq, $config;

    #exit(1);
}
say STDERR "=====> Done.";

exit;
## methods
sub join_pairs {
    my ($sample_size, $forward, $reverse) = @_;

    my $pair_size = sprintf("%.0f", $sample_size/2);
    my $intseq = File::Spec->catfile($opts{output_dir}, $base."_$sample_size".'_interl.fastq.gz');
    my $cmd = "pairfq joinpairs -f <(seqtk sample -s 11 $forward $pair_size) -r <(seqtk sample -s 11 $reverse $pair_size)";
    $cmd .= " -o $intseq -c gzip";

    system_bash($cmd); 

    return $intseq;
}

sub get_mem_usage {
    my $cwd = getcwd();

    my @mems;
    find( sub { push @mems, $File::Find::name if -f and /massif.*/ }, $cwd );

    my $memfile = (nsort @mems)[0];

    open my $in, '-|', 'ms_print', $memfile;

    my ($size, $mem);
    while (my $line = <$in>) {
	chomp $line;
	if ($line =~ /^\s+([MG]B)/ ) {
	    $size = $1;
	    $mem = <$in>;
	    chomp $mem;
	    $mem =~ s/\^.*//;
	}
    }
    close $in;

    unlink @mems;
    return join q{ }, $mem, $size;
}

sub write_results {
    my ($xcov, $sample_size, $logfile, $sumfile, $mem, $out) = @_;

    open my $in, '<', $logfile;

    my ($clstot, $uclstot, $repfrac, $singfrac, $reptot, $annotot, $time);
    while (my $line = <$in>) {
	chomp $line;
	if ($line =~ /Results - Total sequences clustered:\s+($RE{num}{real})/) {
	    $clstot = $1;
	}
	if ($line =~ /Results - Total sequences unclustered:\s+($RE{num}{real})/) {
	    $uclstot = $1;
	}
	if ($line =~ /Results - Repeat fraction from clusters:\s+($RE{num}{real})/) {
	    $repfrac = $1;
	}
	if ($line =~ /Results - Singleton repeat fraction:\s+($RE{num}{real})/) {
	    $singfrac = $1;
	}
	if ($line =~ /Results - Total repeat fraction:\s+($RE{num}{real})/) {
	    $reptot = $1;
	}
	if ($line =~ /Results - Total repeat fraction from annotations:\s+($RE{num}{real})/) {
	    $annotot = $1;
	}
	if ($line =~ /======== Transposome completed at:.*Elapsed time: (.*). ========/) {
	    $time = $1;
	}
    }

    open my $sum, '<', $sumfile;
    
    my $famct = 0;
    while (my $line = <$sum>) {
	chomp $line;
	next if $line =~ /^ReadNum/;
	$famct++;
    }
    close $sum;

    say $out join "\t", "${xcov}X ($sample_size reads)", $clstot, $uclstot, $repfrac, $singfrac, $reptot, $annotot, $famct, $time, $mem;

    return;
}

sub write_config {
    my ($run_opts) = @_;

    my $config = 
"blast_input:
  - sequence_file:      $run_opts->{sequence_file}
  - sequence_format:    $run_opts->{sequence_format}
  - thread:             $run_opts->{threads}
  - output_directory:   $run_opts->{outdir}
clustering_options:
  - in_memory:          1
  - percent_identity:   90
  - fraction_coverage:  0.55
annotation_input:
  - repeat_database:    $run_opts->{repeatdb}
annotation_options:
  - cluster_size:       $run_opts->{cluster_size}
output:
  - run_log_file:       $run_opts->{logfile}
  - cluster_log_file:   $run_opts->{cluster_logfile}";
    
    my $config_file = "transposome_config_$run_opts->{coverage}X.yml";
    open my $out, '>', $config_file;
    say $out $config;
    close $out;

    return $config_file;
}

sub system_bash {
    my @args = ( 'bash', '-c', shift );
    system([0..5], @args);
}

sub usage {
    my $script = basename($0);
    print STDERR <<END

USAGE: $script -s seqs.fastq -o my_cluster_report.txt -n 25000 -c 2 -t 12 [-h]

Required:
 -a|accession            :       A SRA accession. There should be two files found ending in '_1.fastq' and '_2.fastq'.
 -o|output_dir           :       A name for the output to be written.
 -t|threads              :       The number of parallel blast processes to run.
                                 (NB: threads X cpus should be less than or equal to the number
                                 of CPUs on your machine.)
 -repdb|repeat_database  :       A sequence file of repeats to be used for annotation.

Options:
 -f|sequence_format      :       The input sequence format (Default: FASTA).
 -pid|percent_identity   :       Percent identity between pairwise matches in all vs. all blast (Default: 90).
 -fcov|fraction_coverage :       The fraction coverage between pairwise matches in all vs. all blast (Default: 0.55).
 -cls|cluster_size       :       The minimum size of a cluster to be used for annotation (Default: 100).
 -h|help                 :       Print a usage statement.

END
}

