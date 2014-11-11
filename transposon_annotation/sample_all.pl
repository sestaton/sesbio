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
use List::MoreUtils qw(each_array);

my %opt;

GetOptions(\%opt, 'dir|d=s' => \$opt{dir}, 'num|n=i' => \$opt{num}, 'threads|t=i' => \$opt{threads} );

my $cwd = getcwd();

$opt{dir} //= $cwd;
$opt{num} //= 50000;
$opt{threads} //= 1;

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
    #say "DIR: $dir";
    $pm->start($dir) and next;

    my @trimmed_files;
    find( sub { 
	push @trimmed_files, $File::Find::name if -f and /trimmed_screened\.fasta$/; 
          }, $dir );
    
    my @forward = sort grep { /R1/ } @trimmed_files;
    my @reverse = sort grep { /R2/ } @trimmed_files;

    my $ea = each_array(@forward, @reverse);

    while ( my ($f, $r) = $ea->() ) {
	#say join "\n", $f, $r;
	#say "----------------";
	my ($fpfile, $rpfile) = make_pairs($f, $r);
	my $fsample = sample_reads($fpfile, $opt{num});
	my $rsample = sample_reads($rpfile, $opt{num});
	my $intfile = join_pairs($fsample, $rsample);
    }
    $pm->finish;
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
sub make_pairs {
    my ($f, $r) = @_;

    my ($fname, $fpath, $fsuffix) = fileparse($f, qr/\.[^.]*/);
    my ($rname, $rpath, $rsuffix) = fileparse($r, qr/\.[^.]*/);

    my $ffile  = File::Spec->catfile($fpath, $fname."_f".$fsuffix);
    my $rfile  = File::Spec->catfile($rpath, $rname."_r".$rsuffix);
    my @addinfo_f = "pairfq addinfo -i $f -o $ffile -p 1";
    my @addinfo_r = "pairfq addinfo -i $r -o $rfile -p 2";
    run_cmd(\@addinfo_f, 'pairfq addinfo');
    run_cmd(\@addinfo_r, 'pairfq addinfo');
    #say join "\n", $ffile, $rfile;
    #say "------------------------";

    my $fpfile = File::Spec->catfile($fpath, $fname."_fp".$fsuffix);
    my $rpfile = File::Spec->catfile($rpath, $rname."_rp".$rsuffix);
    my $fsfile = File::Spec->catfile($fpath, $fname."_fs".$fsuffix);
    my $rsfile = File::Spec->catfile($rpath, $rname."_rs".$rsuffix);

    my @pairfq_mp = "pairfq makepairs -f $ffile -r $rfile -fp $fpfile -rp $rpfile -fs $fsfile -rs $rsfile";

    run_cmd(\@pairfq_mp, 'pairfq makepairs');
    #say join "\n", $fpfile, $rpfile, $fsfile, $rsfile;
    #say "-----------------------";

    return ($fpfile, $rpfile);
}

sub join_pairs {
    my ($forward, $reverse) = @_;
    
    my ($fname, $fpath, $fsuffix) = fileparse($forward, qr/\.[^.]*/);
    my $name = $fname;
    $name =~ s/R1(?:_)//;
    $name =~ s/fp_//;

    my $ifile  = File::Spec->catfile($fpath, $name."_paired_interl".$fsuffix);
    my @pairfq = "pairfq joinpairs -f $forward -r $reverse -o $ifile";

    run_cmd(\@pairfq, 'pairfq joinpairs');
}

sub sample_reads {
    my ($file, $num) = @_;

    my ($name, $dir, $ext) = fileparse($file, qr/\.[^.]*/);                                                                             
    my $sampled = File::Spec->catfile($dir, $name."_$num".$ext);

    my @seqtk = "seqtk sample -s11 $file $num > $sampled";

    run_cmd(\@seqtk, 'seqtk');

    return $sampled;
}

sub run_cmd {
    my ($cmd, $name) = @_;

    my $exit_value;
    try {
        $exit_value = system([0..5], @$cmd);
    }
    catch {
        say "\nERROR: '$name' exited with exit value: $exit_value. Here is the exception: $_\n";
    };
}
