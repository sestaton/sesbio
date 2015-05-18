#!/usr/bin/env perl

use 5.020;
use strict;
use warnings;
use autodie;
use Cwd;
use File::Find;
use File::Basename;
use Sort::Naturally;
use Parallel::ForkManager;
use IPC::System::Simple qw(system);
use Capture::Tiny       qw(:all);
use Try::Tiny;
use experimental 'signatures';
use Data::Dump;

my @reads;
my $dir      = '/moonriseNFS2/grassa/Highcopy_survey/';
my $samtools = '/home/statonse/github/samtools/samtools';
my $threads  = 1;

find( sub { push @reads, $File::Find::name if -f and /\.rehead.bam$/ }, $dir );

#dd \@reads and exit;

my $pm = Parallel::ForkManager->new($threads);
$pm->run_on_finish( sub { my ($pid, $exit_code, $ident, $exit_signal, $core_dump) = @_;
			  say basename($ident)," just finished with PID $pid and exit code: $exit_code";
		    } );

for my $bam (nsort @reads) {
    my ($acc) = ($bam =~ /(\w+\d+)-?\w+?.rehead\.bam$/);
    $pm->start($acc) and next;
    #say "DEBUG: $acc";
    my ($fq) = bam2fq($bam, $samtools);

    $pm->finish;
    exit;
}
$pm->wait_all_children;

#
# methods
#
sub bam2fq ($bam, $samtools) {
    #say STDERR "===> Converting bam to fastq...";
    
    my $wd = getcwd();
    my ($name, $path, $suffix) = fileparse($bam, qr/\.[^.]*/);
    my $fq  = File::Spec->catfile($wd, $name.".fq");
    my $sfq = File::Spec->catfile($wd, $name."_s.fq");

    my $cmd = "$samtools bam2fq -O -s $sfq $bam > $fq";
    run_cmd($cmd);
    undef $cmd;

    $cmd = "bzip2 $fq";
    run_cmd($cmd);
    undef $cmd;

    $cmd = "bzip2 $sfq";
    run_cmd($cmd);

    return ($fq);
}

sub run_cmd ($cmd) {
    my ($stdout, $stderr, @res) = capture { system([0..5], $cmd); };

    #my @job;
    #try {
	#@job = capture([0..5], $cmd);
    #}
    #catch {
	#say "\nERROR: $cmd exited. Here is the exception: $_\n";
	#exit;
    #};
}

