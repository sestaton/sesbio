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
use IPC::System::Simple qw(system capture);
use Path::Class::File;
use Sort::Naturally;
use Set::IntervalTree;
use Net::SFTP::Foreign;
use Try::Tiny;
use Getopt::Long;
use Data::Dump::Color;
use experimental 'signatures';

my $usage = basename($0).' -o resdir -r reference -u username [-t 8]';

my $host      = 'silo.zoology.ubc.ca';
my $datadir   = File::Spec->catdir('/', 'data', 'raid5part5', 'EST', 'annuus');
my $bwa       = File::Spec->catfile($ENV{HOME}, 'github', 'bwa', 'bwa');
my $samtools  = File::Spec->catfile($ENV{HOME}, 'github', 'samtools', 'samtools');

my %opt;
GetOptions(
    'o|outdir=s'    => \$opt{outdir},
    'u|username=s'  => \$opt{username},
    'r|reference=s' => \$opt{reference},
    't|threads=i'   => \$opt{threads},
    );

die $usage if !$opt{reference} or !$opt{username};

my $elitedir  = 'elite';
my $wilddir   = 'wild';
my $landrcdir = 'landrace';
$opt{threads} //= 8;

#index_ref($bwa, $ref);
my $cwd = getcwd();
my ($refname, $refpath, $refsuf) = fileparse($opt{reference}, qr/\.[^.]*/);
my $sums = $opt{outdir} // File::Spec->catdir($cwd, $refname.'_ests_mapped');

unless (-e $sums) {
    make_path($sums, {verbose => 0, mode => 0711,});
}

for my $dir ($elitedir, $wilddir, $landrcdir) {
    my $sftp = Net::SFTP::Foreign->new($host, user => $opt{username}, autodie => 1);
    $sftp->setcwd($datadir) or die "unable to change cwd: " . $sftp->error;

    say STDERR "LDIR => $dir";
    my ($map) = copy_files($sftp, $dir);
    say STDERR join "map  ====> ";
    dd $map;
    #my $bam = map_reads($opt{reference}, $dir, $bwa, $samtools, $f, $r, $opt{threads});
    #copy $bam, $sums or die "copy failed: $!";
    #remove_tree( $dir, { safe => 1} );
    exit;
}


#
# methods
#
sub map_reads ($ref, $acc, $bwa, $samtools, $f, $r, $threads) {
    my $refbase = $ref;
    $refbase =~ s/\.fa.*//;
    my $prefix = $refbase."_star-aligned";
    my $sbam   = File::Spec->catfile($acc, $acc."_".$refbase."_sorted.bam");

    say STDERR "===> Running STAR on $acc...";
    my $cmd = "$star --runMode alignReads ".
	"--runThreadN 4 ".
	"--genomeDir ".
	"/home/statonse/db ".
	"--readFilesIn $f $r ".
	"--outFileNamePrefix $prefix ".
	"--outSAMtype BAM SortedByCoordinate ".
	"--readFilesCommand zcat ".
	"--outStd BAM_SortedByCoordinate > $bamsort";
    run_cmd($cmd);
    undef $cmd;
    say STDERR "===> Done with STAR.";

    return $sbam;
}

sub copy_files ($sftp, $dir) {
    my $illdir = File::Spec->catdir($dir, 'illumina');
    my $files  = $sftp->ls($illdir,
			   wanted => qr/\.fq$|\.fastq$|\.fastq.gz$|.fq.gz$|.txt$|.txt.gz$/
        ) or die "ls failed: ".$sftp->error;

    unless ( -d $dir ) {
        make_path( $dir, {verbose => 0, mode => 0771,} );
    }
    #dd $files;

    for my $file ( nsort grep { $_->{filename} =~ /_[12].fq$/ } 
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
    
    dd \@fq;
    my (@l, @r, %map);
    my @pairs = nsort @fq;
    my $it    = natatime 2, @pairs;
    while (my @vals = $it->()) {
        my ($fo, $re) = @vals;
	my $base = $fo;
	$base =~ s/_1.f.*//;
        #push @l, $fo;
        #push @r, $re;
	$map{ $base } = join "||", $fo, $re;
    }

    return \%map;
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

sub parse_calls ($calls, $gene_coords) {
    #deletion scf7180038271797:6001-10000 4000 0 1.95573e-05 3.62188e+07 188.22 3.22466e+08 1
    open my $in, '<', $calls;
    while (<$in>) {
	chomp;
	my @f = split;
	my ($ref, $s, $e) = split /[:-]/, $f[1];
	my $res = $gene_coords->{$ref}->fetch($s, $e);
	if ($res) {
	    say join q{ }, @f;
	    dd $res;
	}
    }
}
