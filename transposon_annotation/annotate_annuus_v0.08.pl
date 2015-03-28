#!/usr/bin/env perl

## NB: this is a script for annotating hundreds of sunflower
##     lines with transposome. There are hard-coded paths and
##     queue names, a necessity due to the data storage situation,
##     so this is not going to be generally useful without some
##     modification.
##
## What it do?
##     1) select 2 paired-end files from a breeding line and copy them to a local directory
##     2) downsample the files for efficient screening (seqtk)
##     3) screen out mitochondria, chloroplast and ribosomal repeats (bbsplit.sh)
##     4) downsample the files to 1/2 to total desired coverage for each file (seqtk)
##     5) interleave the files (pairfq)
##     6) run transposome on the cleaned, sampled reads (transposome)
##     7) copy the annotation summary and log to a results directory
##     8) remove all the local results (a single directory) and start another process
##
## TODO: fix the logging method (pass the data to finish -- duh)
use 5.020;
use strict;
use warnings;
use autodie qw(open);
use File::Temp;
use File::Spec;
use File::Find;
use File::Basename;
use File::Copy          qw(copy move);
use File::Path          qw(make_path remove_tree);
use Capture::Tiny       qw(capture_merged);
use IPC::System::Simple qw(system capture);
use Cwd                 qw(getcwd);
use Try::Tiny;
use Parallel::ForkManager;
use experimental 'signatures';

my $data = '/oflow/jmblab/sunflower_1_sequence_data';
my $cwd  = getcwd();
my $sums = File::Spec->catdir($cwd, 'all_lines_summaries_screened_logs');
unless (-e $sums) {
    make_path($sums, {verbose => 0, mode => 0711,});
}

my $outfile = basename($0);
$outfile .= "_all.log";
open my $out, '>>', $outfile or die "\nERROR: Could not open file: $outfile\n";

# these could be options but currently the usage is aimed at being simple
my $threads = 8;
my $matchlen = 50;
my $identity = 50;
my $first_sample = 7e6;
my $sample_level = 1e6;

my @dirs;
find( sub {
    push @dirs, $_ if -d and /^PPN\d+/
    }, $data);

my $pm = Parallel::ForkManager->new($threads);
$pm->run_on_finish( sub { my ($pid, $exit_code, $ident, $exit_signal, $core_dump) = @_;
			  say $out basename($ident)," just finished with PID $pid and exit code: $exit_code";
		      } );

for my $localdir (sort @dirs) {
    $pm->start($localdir) and next;
    my $dirpath = File::Spec->catdir($data, $localdir);
    next unless -e $dirpath;
    my ($line) = ($localdir =~ /PPN(\d+)/);

    my @seqs;
    find( sub {
	push @seqs, $_ if -f and /\.fq$|\.fastq$|\.fastq.gz$|.fq.gz$/
	}, $dirpath);
    
    next unless @seqs;
    my ($f, $r) = (sort @seqs)[0..1];
    my ($fsamp, $rsamp) = sample_seqs($dirpath, $localdir, $f, $r, $first_sample, 1);
    my ($fscr, $rscr) = screen_reads($fsamp, $rsamp, $localdir, $matchlen, $identity);
    my ($fsamp_scr, $rsamp_scr) = sample_seqs($dirpath, $localdir, $fscr, $rscr, $sample_level, 0);
    my $joined = join_reads($localdir, $fsamp_scr, $rsamp_scr, $sample_level);
    my ($summary, $log) = run_transposome($sums, $localdir, $joined, $sample_level);

    $pm->finish;
}

$pm->wait_all_children;
close $out;

#
# methods
#
sub run_transposome ($sums, $localdir, $joined, $readnum) {
    chdir $localdir;
    
    my ($config, $outdir) = write_transposome_conf($localdir, $joined, $readnum);
    unless ( -e $outdir ) {
	make_path( $outdir, {verbose => 0, mode => 0771,} );
    }

    my $script = make_script('transposome');
    open my $tout, '>', $script;
    say $tout "#!/bin/bash\ntransposome --config $config";
    close $tout;
        
    my @job;
    try {
	@job = capture([0..5], "qsub -q rcc-m128-30d -l mem_total=200g $script");
    }
    catch {
	say "\nERROR: transposome exited. Here is the exception: $_\n";
	exit;
    };
    
    my $done = check_job(\@job, $script);
    my @out = glob "$script*";
    unlink @out;
    unlink $config;
    
    my (@sum, @log);
    find( sub {
	push @sum, $File::Find::name if -f and /summary.tsv$/
	}, $outdir);

    @sum = sort @sum;
    my $sumfile = shift @sum;
    move $sumfile, $sums or die "move failed: $!";

    find( sub {
	push @log, $File::Find::name if -f and /log.*\.txt$/
        }, $outdir);

    my $logfile = shift @log;
    move $logfile, $sums or die "move failed: $!";
    remove_tree($outdir);
    chdir "..";
    
    # remove the data from this line
    remove_tree($localdir);
    return ($sumfile, $logfile);
}

sub check_job ($job, $script) {
    my $id; 
    for (@$job) {
	if (/(\d+)/) {
	    $id = $1;
	}
    }

    die "\nERROR: $job can not be found. Check $script."
	unless defined $id;
    my $scrname = basename($script);
    my $cmd = "qstat -j $id";

    attempt : {
	my $out = capture_merged { system([0..5], $cmd); };
	
	if ($out =~ /jobs do not exist/) {
	    return 1;
	}
	else {
	    sleep 60;
	    redo attempt;
	}
    }
}

sub write_transposome_conf ($dir, $joined, $readnum) {
    my $config = "transposome_config_".$dir.".yml";
    my $outdir = $dir."_transposome_results_out";

    open my $out, '>', $config;
    say $out "
blast_input:
  - sequence_file:      $joined
  - sequence_format:    fastq
  - sequence_num:       50000
  - cpu:                1
  - thread:             1
  - output_directory:   $outdir
clustering_options:
  - in_memory:          1
  - percent_identity:   90
  - fraction_coverage:  0.55
  - merge_threshold:    1000
annotation_input:
  - repeat_database:    ~/db/RepBase1801_Hannuus_LTR-RT_families_repbase_frmt.fasta
annotation_options:
  - cluster_size:       100
  - blast_evalue:       10
output:
  - run_log_file:       ${dir}_${readnum}_log.txt
  - cluster_log_file:   ${dir}_${readnum}_cluster_report.txt";

    close $out;

    return ($config, $outdir);
}

sub join_reads ($dir, $fsamp, $rsamp, $sample_level) {
    my $joined = $dir."_${sample_level}_interl.fastq";
    my $jout   = File::Spec->catfile($dir, $joined);

    if (defined $fsamp && defined $rsamp && -s $fsamp && -s $rsamp) {
        my $jcmd = "pairfq joinpairs -f $fsamp -r $rsamp -o $jout";

        my $script = make_script('join');
        open my $pjout, '>', $script;
        say $pjout "#!/bin/bash\n$jcmd";
        close $pjout;

	my @job;
	try {
	    @job = capture([0..5], "qsub -q rcc-30d $script");
	}
	catch {
	    say "\nERROR: pairfq exited. Here is the exception: $_\n";
	    exit;
	};
	
	my $done = check_job(\@job, $script);
        my @out = glob "$script*";
        unlink @out, $fsamp, $rsamp;
	
	return $joined;
    }
}

sub sample_seqs ($dirpath, $dir, $f, $r, $level, $addinfo) {
    unless (-e $dir) {
        make_path($dir, {verbose => 0, mode => 0711,});
    }

    if ($addinfo) {
	my $ext = ".fastq";
	my $pair_level = $level / 2;
	my $forig   = File::Spec->catfile($dirpath, $f);
	my $rorig   = File::Spec->catfile($dirpath, $r);
	my $forward = File::Spec->catfile($dir, $f);
	my $reverse = File::Spec->catfile($dir, $r);
	
	my ($fpair, $rpair) = ($forward, $reverse);
	$fpair =~ s/\.gz// if $fpair =~ /\.gz$/;
	$rpair =~ s/\.gz// if $rpair =~ /\.gz$/;
	my ($ffile, $fdir, $fext) = fileparse($fpair, qr/\.[^.]*/);
	my ($rfile, $rdir, $rext) = fileparse($rpair, qr/\.[^.]*/);
	
	$fpair    = File::Spec->catfile($dir, $ffile."_$pair_level".$ext);
	$rpair    = File::Spec->catfile($dir, $rfile."_$pair_level".$ext);
	my $fsamp = File::Spec->catfile($dir, $ffile."_1_$pair_level".$ext);
	my $rsamp = File::Spec->catfile($dir, $rfile."_2_$pair_level".$ext);
	
	my ($script, $done, @job, @out);

	$script = make_script('copy');
	open my $out, '>', $script;
	say $out "#!/bin/bash\ncp $forig $forward\ncp $rorig $reverse";
	close $out;
	
	try {
	    @job = capture([0..5], "qsub -q copyq $script");
	}
	catch {
	    say "\nERROR: Copy failed: $_\n";
	};
	
	$done = check_job(\@job, $script);
	@out = glob "$script*";
	unlink @out;
	undef $script;
	undef @out;
	undef $done;
	undef @job;

	my $fsample_cmd = "zcat $forward | seqtk sample -s11 - $pair_level | pairfq addinfo -i - -o $fsamp -p 1";
	my $rsample_cmd = "zcat $reverse | seqtk sample -s11 - $pair_level | pairfq addinfo -i - -o $rsamp -p 2";

	$script = make_script('sample');
	open my $saout, '>', $script;
	say $saout "#!/bin/bash\n$fsample_cmd\n$rsample_cmd";
	close $saout;
	
	try {
	    @job = capture([0..5], "qsub -q rcc-30d $script");
	}
	catch {
	    say "\nERROR: seqtk exited. Here is the exception: $_\n";
	    exit;
	};
	
	$done = check_job(\@job, $script);
	@out = glob "$script*";
	unlink @out;
	undef @out;
	undef @job;
    
	unlink @out, $forward, $reverse, $fpair, $rpair;
	
	return ($fsamp, $rsamp);
    }
    else {
	my $ext = ".fastq";
        my $pair_level = $level / 2;
	my ($forward, $reverse) = ($f, $r);
        my ($ffile, $fdir, $fext) = fileparse($forward, qr/\.[^.]*/);
        my ($rfile, $rdir, $rext) = fileparse($reverse, qr/\.[^.]*/);

        my $fpair = File::Spec->catfile($dir, $ffile."_$pair_level".$ext);
        my $rpair = File::Spec->catfile($dir, $rfile."_$pair_level".$ext);

	my $fsample_cmd = "seqtk sample -s11 $forward $pair_level > $fpair";
        my $rsample_cmd = "seqtk sample -s11 $reverse $pair_level > $rpair";

        my $script = make_script('sample');
	open my $saout, '>', $script;
        say $saout "#!/bin/bash\n$fsample_cmd\n$rsample_cmd";
        close $saout;

	my @job;
        try {
            @job = capture([0..5], "qsub -q rcc-30d $script");
        }
        catch {
            say "\nERROR: seqtk exited. Here is the exception: $_\n";
            exit;
        };

        my $done = check_job(\@job, $script);
        my @out = glob "$script*";
	unlink @out;

	unlink $forward, $reverse;
	return ($fpair, $rpair);
    }
}

sub screen_reads ($fsamp, $rsamp, $dir, $matchlen, $identity) {
    my ($ffile, $fdir, $fext) = fileparse($fsamp, qr/\.[^.]*/);
    my ($rfile, $rdir, $rext) = fileparse($rsamp, qr/\.[^.]*/);
    my $contam = File::Spec->catfile($dir, $ffile."_contam_%".".fastq");
    my $fsamp_scrfq = File::Spec->catfile($dir, $ffile."_scr".".fastq");
    my $rsamp_scrfq = File::Spec->catfile($dir, $rfile."_scr".".fastq");

    ($fsamp_scrfq, $rsamp_scrfq) = run_bbsplit($dir, $fsamp, $rsamp, $contam, $fsamp_scrfq, $rsamp_scrfq);

    return ($fsamp_scrfq, $rsamp_scrfq);
}

sub run_bbsplit ($dir, $fsamp, $rsamp, $contam, $fsamp_scrfq, $rsamp_scrfq) {
    my $db = '/home/jmblab/statonse/db/ScreenDB_Hannus-cpDNA_UniVec-5-2_Plant-mtDNA_LSU-SSU-DB.fasta';
    my $bbcmd = "~/apps/bbmap/bbsplit.sh -Xmx8g in1=$fsamp in2=$rsamp ref=$db ";
    $bbcmd .= "basename=$contam outu1=$fsamp_scrfq outu2=$rsamp_scrfq";

    my $bbscript = make_script('bbsplit');

    open my $out, '>', $bbscript;
    say $out "#!/bin/bash\n$bbcmd";
    close $out;

    my @job;
    try {
	@job = capture([0..5], "qsub -q rcc-30d -l mem_total=8g $bbscript");
    }
    catch {
	say "\nERROR: blast exited. Here is the exception: $_\n";
	exit;
    };

    my $done = check_job(\@job, $bbscript);
    undef $done;
    my @bout = glob "$bbscript*";
    my @out = glob "$dir/*ScreenDB*";
    unlink @bout, @out;

    return ($fsamp_scrfq, $rsamp_scrfq);
}

sub make_script ($command) {
    my $cwd = getcwd();
    my $tmpiname = $command."_XXXX";
    my $fname = File::Temp->new( TEMPLATE => $tmpiname,
				 DIR      => $cwd,
				 SUFFIX   => ".sh",
				 UNLINK   => 0);
    
    my $script = $fname->filename;

    return $script;
}
