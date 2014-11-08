#!/usr/bin/env perl

use 5.020;
use strict;
use warnings;
use File::Find;
use File::Spec;
use File::Basename;
use IPC::System::Simple qw(system);
use Try::Tiny;
use Getopt::Long;
use POSIX qw(strftime);
use Time::HiRes qw(gettimeofday);
use Parallel::ForkManager;
use Cwd;

my %opt;

GetOptions(\%opt, 'dir|d=s' => \$opt{dir}, 'num|n=i' => \$opt{num}, 'threads|t=i' => \$opt{threads} );

my $cwd = getcwd();

$opt{dir} //= $cwd;
$opt{num} //= 50000;
$opt{threads} //= 12;

my $t1 = gettimeofday();
my $fs = POSIX::strftime('%d-%m-%Y %H:%M:%S', localtime);
say "***** Begin sampling reads at: $fs.";

my $pm = Parallel::ForkManager->new( $opt{threads} );
$pm->run_on_finish( sub { my ($pid, $exit_code, $ident, $exit_signal, $core_dump) = @_;
			  my $t2 = gettimeofday();
			  my $elapsed = $t2 - $t1;
			  my $time = sprintf("%.2f",$elapsed/60);
			  say basename($ident),
			  " just finished with PID $pid and exit code: $exit_code in $time minutes";
		      } );

my @dirs;
find( sub { 
    push @dirs, $File::Find::name if -d and /sample/i or /moscow/i; 
      }, $opt{dir} );

for my $dir (@dirs) {
    my @trimmed_files;
    find( sub { 
	push @trimmed_files, $File::Find::name if -f and /trimmed_screened\.fasta$/; 
          }, $dir );

    for my $file (@trimmed_files) {
	my ($name, $dir, $ext) = fileparse($file, qr/\.[^.]*/);
	my $sampled = File::Spec->catfile($dir, $name."_50k".$ext);
	$pm->start($sampled) and next;
	#warn "File exists: $sampled, skipping..." and next if -e $sampled;
	sample_reads($file, $sampled, $opt{num});
	$pm->finish;
    }
}

$pm->wait_all_children;

my $t2 = gettimeofday();
my $total_elapsed = $t2 - $t1;
my $final_time = sprintf("%.2f",$total_elapsed/60);
my $ft = POSIX::strftime('%d-%m-%Y %H:%M:%S', localtime);
say "***** Finished sampling reads in $final_time minutes at: $ft.";

#
# methods
#
sub sample_reads {
    my ($file, $sampled, $num) = @_;

    my @seqtk = "seqtk sample -s11 $file $num > $sampled";

    my $exit_value;
    try {
	$exit_value = system([0..5], @seqtk);
    }
    catch {
	say "\nERROR: 'seqtk' exited with exit value: $exit_value. Here is the exception: $_\n";
    };
}
