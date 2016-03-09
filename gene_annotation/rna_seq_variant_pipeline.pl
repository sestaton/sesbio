#!/usr/bin/env perl

use 5.022;
use strict;
use warnings;
use Cwd;
use File::Find;
use File::Spec;
use File::Basename;
use File::Path          qw(make_path remove_tree);
use File::Copy          qw(copy);
use List::MoreUtils     qw(natatime);
use IPC::System::Simple qw(system);
use Capture::Tiny       qw(capture);
use Path::Class::File;
use Sort::Naturally;
use Set::IntervalTree;
use Net::SFTP::Foreign;
use Parallel::ForkManager;
use Try::Tiny;
use Getopt::Long;
use Data::Dump::Color;
use experimental 'signatures';

my $usage = basename($0).' -h host.name -d datadir -o resdir -r reference -u username [-t 8]';

my $star     = File::Spec->catfile($ENV{HOME}, 'github', 'STAR', 'bin', 'Linux_x86_64', 'STAR');
my $samtools = File::Spec->catfile($ENV{HOME}, 'github', 'samtools', 'samtools');
my $java     = File::Spec->catfile('/', 'etc', 'alternatives', 'jre_1.8.0', 'bin', 'java');
my $picard   = File::Spec->catfile($ENV{HOME}, 'apps', 'picard-tools-2.1.0', 'picard.jar');
my $gatk     = File::Spec->catfile($ENV{HOME}, 'apps', 'bin', 'GenomeAnalysisTK.jar');
my $frbayes  = File::Spec->catfile($ENV{HOME}, 'github', 'freebayes', 'bin', 'freebayes');

my %opt;
GetOptions(
    'h|host=s'      => \$opt{host},
    'd|datadir=s'   => \$opt{datadir},
    'o|outdir=s'    => \$opt{outdir},
    'u|username=s'  => \$opt{username},
    'r|reference=s' => \$opt{reference},
    't|threads=i'   => \$opt{threads},
    );

die $usage if !$opt{reference} or !$opt{username} or !$opt{host} or !$opt{datadir};

my $elitedir  = 'elite';
my $wilddir   = 'wild';
my $landrcdir = 'landrace';
$opt{threads} //= 2;

my $cwd = getcwd();
my ($refname, $refpath, $refsuf) = fileparse($opt{reference}, qr/\.[^.]*/);
my $sums;
if ($opt{outdir}) {
    $sums = File::Spec->catdir($cwd, $opt{outdir});
}
else {
    $sums = File::Spec->catdir($cwd, $refname.'_gatk_rnaseq_vars');
}

unless (-e $sums) {
    make_path($sums, {verbose => 0, mode => 0711,});
}

for my $dir ($elitedir, $wilddir, $landrcdir) {
    my $dir_sums = File::Spec->catdir($sums, $dir);
    unless (-e $dir_sums) {
        make_path($dir_sums, {verbose => 0, mode => 0711,});
    }

    say STDERR "===> Transferring data for: $dir";
    my $sftp = Net::SFTP::Foreign->new($opt{host}, user => $opt{username}, autodie => 1);
    $sftp->setcwd($opt{datadir}) or die "unable to change cwd: " . $sftp->error;
    my ($map) = copy_files($sftp, $dir);
    
    say STDERR "===> Mapping results for: $dir";
    map_reads($refname, $opt{reference}, $map, $dir, $dir_sums,
	      $star, $samtools, $java, $picard, $gatk, $opt{threads});
    chdir $cwd or die $!;
    remove_tree( $dir, { safe => 1} ); ## remove sequence data
}

chdir $cwd or die $!;
my $merged_bams = merge_bams($java, $picard, $sums, $elitedir, $wilddir, $landrcdir);
my $vcfs = call_variants($frbayes, $opt{reference}, $merged_bams);
#
# methods
#
sub merge_bams ($java, $picard, $sums, $elitedir, $wilddir, $landrcdir) {
    my %bams;
    for my $dir ($elitedir) {
	my $dir_sums = File::Spec->catdir($sums, $dir);
	my @files;
	find( sub { push @files, $File::Find::name if -f and /\.bam$/ }, $dir_sums);
	if (@files > 0) {
	    my $inc_str;
	    for my $f (@files) {
		$inc_str .= "I=$f ";
	    }
	    my $merged_bam = File::Spec->catfile($sums, $dir, $dir.'_merged.bam');
	    
	    my $cmd = "$java -jar $picard MergeSamFiles ".
		"CREATE_INDEX=true ".
		"$inc_str ".
		"O=$merged_bam";
	    run_cmd($cmd);
	    $bams{ $dir } = $merged_bam;
	}
	else {
	    say STDERR "\nERROR: No BAM files found in '$dir'.\n";
	}
    }
    return \%bams;
}

sub call_variants ($frbayes, $ref, $bams) {
    my %vcfs;
    for my $acc (keys %$bams) {
	my $bam = $bams->{$acc};
	my ($file, $path, $ext) = fileparse($bam, qr/\.[^.]*/);
	my $vcf = File::Spec->catfile($path, $file.'_vars.vcf');

	my $cmd = "$frbayes -f $ref $bam > $vcf";
	run_cmd($cmd);
	$vcfs { $acc }  = $vcf;
    }
    return \%vcfs;
}

sub map_reads ($refname, $ref, $map, $dir, $dir_sums, $star, $samtools, $java, $picard, $gatk, $threads) {
    chdir $dir or die $!;
    my $wd = getcwd();

    my $pm = Parallel::ForkManager->new($threads);
    local $SIG{INT} = sub {
        warn "Caught SIGINT; Waiting for child processes to finish.";
        $pm->wait_all_children;
        exit 1;
    };

    $pm->run_on_finish( sub { 
	my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data_structure_reference) = @_;
	if (defined $data_structure_reference) {
	    my ($acc, $bamfile) = @$data_structure_reference;
	    say STDERR "===> Finished processing $acc (PID: $pid) with exit signal: $exit_signal.";
	    copy $bamfile, $dir_sums or die $!;
	    remove_tree( $acc, { safe => 1} );
	}
	else { 
	    say STDERR "WARN: No message received from child process $pid!\n";
	}
      } 
    );

    my (%bams, %vcfs);
    for my $acc (nsort keys %$map) {
	next if $acc =~ /ha412/i;
	$pm->start($acc) and next;
	$SIG{INT} = sub { $pm->finish };

	my $bamfile = align_and_filter($ref, $refname, $map, $acc, $wd, 
				       $star, $samtools, $java, $picard, $gatk);
	my @res = ($acc, $bamfile);
	$pm->finish(0, \@res);
    }
    $pm->wait_all_children;

    return;
}

sub align_and_filter ($ref, $refname, $map, $acc, $wd, $star, $samtools, $java, $picard, $gatk) {
    unless (-e $acc) {
	make_path($acc, {verbose => 0, mode => 0711,});
    }

    chdir $acc or die $!;
    my ($f, $r) = split /\|\|/, $map->{$acc};
    my $prefix = $refname."_".$acc."_star-aligned";
    my $sam    = $acc."_".$refname.".sam";
    my $sam2   = $acc."_".$refname."_pass2.sam";
    my $stardb = $acc."_stargenome";
    my $sjtab  = $prefix."SJ.out.tab";
    my $dbdir  = File::Spec->catdir($ENV{HOME}, 'db'); # location of STAR index

    unless (-e $stardb) {
	make_path($stardb, {verbose => 0, mode => 0711,});
    }

    ## STAR 1st pass
    say STDERR "===> Running STAR on $acc...";
    my $starcmd = "$star --runMode alignReads ".
	"--runThreadN 4 ".
	"--genomeDir $dbdir ".
	"--readFilesIn ../$f ../$r ".
	"--outFileNamePrefix $prefix ".
	"--readFilesCommand zcat ".
	"--outStd SAM > $sam";
    my @star1 = capture_cmd($starcmd);
    undef $starcmd;
    
    ## STAR 2nd pass # todo: redirect output
    $starcmd = "$star --runMode genomeGenerate ".
	"--genomeDir $stardb ".
	"--genomeFastaFiles $ref ".
	"--sjdbFileChrStartEnd $sjtab ".
	"--sjdbOverhang 75 ".
	"--runThreadN 4";
    run_cmd($starcmd);
    undef $starcmd;
    
    $starcmd = "$star --runMode alignReads ".
	"--genomeDir $stardb ".
	"--readFilesIn ../$f ../$r ".
	"--readFilesCommand zcat ".
	"--runThreadN 4 ".
	"--outStd SAM > $sam2";
    run_cmd($starcmd);
    undef $starcmd;

    ## filtering
    my $rg_bam  = File::Spec->catfile($wd, $acc, $acc."_rg_added_sorted.bam");
    my $rg_dbam = File::Spec->catfile($wd, $acc, $acc."_rg_added_sorted_dedupped.bam"); 
    my $rg_sbam = File::Spec->catfile($wd, $acc, $acc."_rg_added_sorted_dedupped_split.bam");
    my $rg_log  = File::Spec->catfile($wd, $acc, $acc."_cigarsplitn.log");

    my $piccmd = "$java -Xmx10g -jar $picard AddOrReplaceReadGroups ".
	"I=$sam2 ".
	"O=$rg_bam ".
	"SO=coordinate RGID=id_$acc RGLB=lib_$acc RGPL=illumina RGPU=machine_$acc RGSM=$acc"; 
    my ($astdo, $astderr, @ars) = capture_cmd($piccmd);
    undef $piccmd;

    $piccmd = "$java -Xmx10g -jar $picard MarkDuplicates ".
	    "I=$rg_bam ".
	    "O=$rg_dbam ".
	    "CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics";
    my @md = capture_cmd($piccmd);
    undef $piccmd;

    $piccmd = "$java -Xmx10g -jar $gatk ".
	    "-T SplitNCigarReads ".
	    "-R $ref ".
	    "-I $rg_dbam ".
	    "-o $rg_sbam ".
	    "-rf ReassignOneMappingQuality ".
	    "-RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS";
    my $err = 0;
    my ($stdout, $stderr, $out) = capture_cmd($piccmd);
    open my $el, '>>', $rg_log;
    say $el $stdout;
    say $el $stderr;
    close $el;

    for my $l (split /^/, $stderr) {
	chomp $l;
	if ($l =~ /error.*?encoding.*?high/i) {
	    $piccmd .= " --fix_misencoded_quality_scores";
	    $err = 1;
	}
    }
    my @spln = capture_cmd($piccmd) if $err;
    undef $piccmd;

    ## indel realignment
    my $ind_int = File::Spec->catfile($wd, $acc, $acc."_indel.intervals");
    my $indcd = "$java -Xmx15g -jar $gatk ".
	"-T RealignerTargetCreator ".
	"-R $ref ".
	"-o $ind_int ".
	"-I $rg_sbam";
    my @realn = capture_cmd($indcd);
    undef $indcd;

    my $ira_bam = File::Spec->catfile($wd, $acc, $acc."_indelRealigned.bam");
    my $racmd = "$java -Xmx10g -jar $gatk ".
	"-I $rg_sbam ".
	"-R $ref ".
	"-T IndelRealigner ".
	"-targetIntervals $ind_int ".
	"-o $ira_bam";
    my @indrealn = capture_cmd($racmd);
    undef $racmd;

    return $ira_bam;
}

sub copy_files ($sftp, $dir) {
    my $illdir = File::Spec->catdir($dir, 'illumina');
    my $files  = $sftp->ls($illdir,
			   wanted => qr/\.fq$|\.fastq$|\.fastq.gz$|.fq.gz$|.txt$|.txt.gz$/
        ) or die "ls failed: ".$sftp->error;

    unless ( -d $dir ) {
        make_path( $dir, {verbose => 0, mode => 0771,} );
    }

    for my $file ( nsort grep { $_->{filename} =~ /_[12]\.(?:fq|txt)\.gz$/ } 
		   grep { $_->{longname} !~ /^l/ } @$files ) {
        my $flocal = File::Spec->catfile($dir, $file->{filename});
	say STDERR "FILE => $flocal";
        my $size = $file->{a}->size;
	
        my $rem_f = File::Spec->catfile($datadir, $illdir, $file->{filename});
	say STDERR "REMFILE => $rem_f";
        $sftp->get($rem_f, $flocal) or die "get failed: " . $sftp->error;
        
        my $lfsize = -s $flocal;
        die "Failed to fetch complete file: $flocal (local size: $lfsize, remote size: $size)"
            unless $size == $lfsize;
    }
    
    my @fq;
    find( sub { push @fq, $File::Find::name 
                    if -f and /\.fq$|\.fastq$|\.fastq.gz$|.fq.gz$|.txt$|.txt.gz$/ }, $dir );
    
    my (@l, @r, %map);
    my @pairs = nsort @fq;
    my $it    = natatime 2, @pairs;
    while (my @vals = $it->()) {
        my ($fo, $re) = @vals;
	my $base = $fo;
	$base =~ s/_1.(?:f|t).*//;
	$base =~ s/$dir\///;
	$fo   =~ s/$dir\///;
	$re   =~ s/$dir\///;

	$map{ $base } = join "||", $fo, $re;
    }

    return \%map;
}

sub capture_cmd ($cmd) {
    my ($stdout, $stderr, @out) = capture { system([0..5], $cmd) };

    return ($stdout, $stderr, \@out);
}

sub run_cmd ($cmd) {
    my @job;
    try {
        @job = system([0..5], $cmd);
    }
    catch {
        say "\nERROR: $cmd exited. Here is the exception: $_\n";
    };
}
