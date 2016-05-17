#!/usr/bin/env perl

use 5.022;
use strict;
use warnings;
use autodie;
use File::Find;
use File::Spec;
use File::Basename;
use File::Path          qw(make_path);
use IPC::System::Simple qw(system);
use Set::IntervalTree;
use Try::Tiny;
use Getopt::Long;
use Data::Dump::Color;
use experimental 'signatures';

my $usage    = basename($0).' -i bamdir -o cnvresdir -f ref_split_fas_dir';
my $cnvnator = join "/", $ENV{HOME}, 'apps', 'CNVnator_v0.3.2', 'src', 'cnvnator';
my $chrom    = 'scf7180038271797 scaffold_500';

my %opt;
GetOptions(
    'i|indir=s'    => \$opt{indir},
    'o|outdir=s'   => \$opt{outdir},
    'f|fastadir=s' => \$opt{fasdir},
    'b|binsize=i'  => \$opt{binsize},
    );

die $usage if !$opt{indir} or !$opt{outdir} or !$opt{fasdir};

unless ( -d $opt{outdir} ) {
    make_path( $opt{outdir}, {verbose => 0, mode => 0771,} );
}

$opt{binsize} //= 100;

my $scaf500_tree = Set::IntervalTree->new;
my $scf7180_tree = Set::IntervalTree->new;
$scaf500_tree->insert('lg4', 29100, 36600);
$scf7180_tree->insert('lg16', 33300, 38600);

my %gene_coords = (
    'scf7180038271797' => $scf7180_tree,
    'scaffold_500'     => $scaf500_tree,
    );

my $rootfile = $opt{indir}."_cnvnator_merged.root";
my $callfile = $opt{indir}."_cnvnator_merged.calls";

my (@bams, @roots);
find( sub { push @bams, $File::Find::name if -f and /\.bam$/ }, $opt{indir} );
my $bams = join q{ }, @bams;
#say $bams;
for my $bam (@bams) {
    my ($name, $path, $suffix) = fileparse($bam, qr/\.[^.]*/);
    my $rootfile = File::Spec->catfile($opt{outdir}, $name.".root");
    my $mkroot_cmd = "$cnvnator -root $rootfile -chrom $chrom -tree $bams";
    run_cmd($mkroot_cmd);
}

find( sub { push @roots, $File::Find::name if -f and /\.root$/ }, $opt{outdir} );
my $roots = join q{ }, @roots;
my $mkroot_cmd = "$cnvnator -root test_merged.root -chrom $chrom -merge $roots";
run_cmd($mkroot_cmd);
my $makehist_cmd = "$cnvnator -root test_merged.root -chrom $chrom -his $opt{binsize} -d $opt{fasdir}";
run_cmd($makehist_cmd);
my $stat_cmd = "$cnvnator -root test_merged.root -chrom $chrom -stat $opt{binsize} -d $opt{fasdir}";
run_cmd($stat_cmd);
my $part_cmd = "$cnvnator -root test_merged.root -chrom $chrom -partition $opt{binsize} -d $opt{fasdir}";
run_cmd($part_cmd);
my $call_cmd = "$cnvnator -root test_merged.root -chrom $chrom -call $opt{binsize} -d $opt{fasdir} > test_merged.calls";
run_cmd($call_cmd);
parse_calls($callfile, \%gene_coords);

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
	if (@$res) {
	    say "CNV on target region: ",@$res;
	    say join q{ }, @f;
	}
    }
}
