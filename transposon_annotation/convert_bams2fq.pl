#!/usr/bin/env perl

use 5.020;
use warnings;
use autodie;
use Cwd;
use File::Find;
use File::Basename;
use Sort::Naturally;
use Parallel::ForkManager;
use IPC::System::Simple qw(system);
use Capture::Tiny       qw(:all);
use experimental 'signatures';
#use Data::Dump::Color;

my @reads;
my $usage = "perl ".basename($0)." dir_with_bams\n";
my $dir   = shift or die $usage;      
my $samtools = findprog('samtools');
my $threads  = 1;

find( sub { push @reads, $File::Find::name if -f and /\.rehead.bam\z/ }, $dir );

my $pm = Parallel::ForkManager->new($threads);
$pm->run_on_finish( sub { my ($pid, $exit_code, $ident, $exit_signal, $core_dump) = @_;
			  say basename($ident)," just finished with PID $pid and exit code: $exit_code";
		    } );

for my $bam (nsort @reads) {
    my ($acc) = ($bam =~ /(\w+\d+)-?\w+?.rehead\.bam\z/);
    $pm->start($acc) and next;
    my ($fq) = bam2fq($bam, $samtools);

    $pm->finish;
    exit;
}
$pm->wait_all_children;

#
# methods
#
sub bam2fq ($bam, $samtools) {
    my $wd = getcwd();
    my ($name, $path, $suffix) = fileparse($bam, qr/\.[^.]*/);
    my $fq  = File::Spec->catfile($wd, $name.".fq.bz2");
    my $sfq = File::Spec->catfile($wd, $name."_s.fq");

    my $cmd = "$samtools bam2fq -O -s $sfq $bam | bzip2 > $fq";
    run_cmd($cmd);
    undef $cmd;

    $cmd = "bzip2 $sfq";
    run_cmd($cmd);

    return ($fq);
}

sub run_cmd ($cmd) {
    my ($stdout, $stderr, @res) = capture { system([0..5], $cmd); };
}

sub findprog ($prog) {
    my ($path, $err, @res) = capture { system([0..5], "which $prog 2> /dev/null") };    
    chomp $path;
    if ( ! -e $path && ! -x $path ) {
        die "\nERROR: Cannot find $prog binary. Exiting.\n\n";
    } 
    else {
        return $path;
    }
}
