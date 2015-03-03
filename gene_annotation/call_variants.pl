#!/usr/bin/env perl

use 5.020;
use strict;
use warnings;
use autodie;
use File::Find;
use File::Basename;
use Parallel::ForkManager;
use IPC::System::Simple qw(system);
use Capture::Tiny       qw(:all);
use Try::Tiny;
use List::MoreUtils     qw(natatime);
use experimental 'signatures';
use Data::Dump;

my @reads;
my $usage    = "$0 ref\n";
my $ref      = shift or die $usage;
my $species  = 'scf7180038271797';
my $dir      = 'cleaned_transcriptome_reads';
my $bwa      = '/home/statonse/github/bwa/bwa';
my $samtools = '/home/statonse/github/samtools/samtools';
my $bcftools = '/home/statonse/github/bcftools/bcftools';
my $vcfutils = '/home/statonse/github/bcftools/vcfutils.pl';
my $threads  = 2;

find( sub { push @reads, $File::Find::name if -f and /\.fq$/ }, $dir );
my @rawreads = grep { /[12]/ } 
               grep { ! /EPSPS_reads|unpe|test/ } @reads;

#dd \@rawreads and exit;

my @sreads = sort @rawreads;
my $it = natatime 2, @sreads;

my $pm = Parallel::ForkManager->new(12);
$pm->run_on_finish( sub { my ($pid, $exit_code, $ident, $exit_signal, $core_dump) = @_;#
			  say basename($ident)," just finished with PID $pid and exit code: $exit_code";
		    } );

while (my @vals = $it->()) {
    my ($f, $r) = @vals;
    my ($acc) = ($f =~ /[weedy|wild|elite]\.(.*)\.[12]\.fq$/);
    #weedy.Academy2.1.fq
    #@say "Debug: $f $r" and exit unless defined $acc;
    $pm->start($acc) and next;
    #say $acc;

    #index_ref($ref, $bwa, $samtools);
    #say "DEBUG: $acc";
    my $bam = map_reads($ref, $acc, $bwa, $f, $r, $threads);
    my ($fvcf, $cvcf) = call_variants($ref, $bam, $samtools, $bcftools, $vcfutils);
    run_vep($fvcf, $cvcf, $species);

    $pm->finish;
    #exit;
}
$pm->wait_all_children;

#
# methods
#
sub call_variants ($ref, $bam, $samtools, $bcftools, $vcfutils) {
    #say STDERR "===> Calling variants...";
    
    my $bcf = $bam;
    $bcf =~ s/\.bam//;
    my $vcf = $bcf;
    my $calls = $bcf;
    $bcf .= "_raw.bcf";
    $vcf .= "_filt.vcf";
    $calls .= "_calls.vcf";

    my $cmd = "$samtools mpileup -d8000 -uf $ref $bam | $bcftools ";
    $cmd .= "view -c 0 -g \"^miss\" - > $bcf";
    run_cmd($cmd);
    undef $cmd;

    $cmd = "$bcftools view $bcf | $vcfutils varFilter -D100 > $vcf";
    run_cmd($cmd);
    undef $cmd;

    # calls
    $cmd = "$samtools mpileup -d8000 -uf $ref $bam | $bcftools call -mv > $calls";
    run_cmd($cmd);

    return ($vcf, $calls);
}

sub map_reads ($ref, $acc, $bwa, $f, $r, $threads) {
              #$ref, $acc, $bwa, $f, $r, $threads);
    #say STDERR "===> Mapping reads $f and $r to $ref ...";
    my $refbase = $ref;
    $refbase =~ s/\.fa.*//;
    my $sam  = $acc."_".$refbase.".sam";
    my $bam  = $acc."_".$refbase.".bam";
    my $sort = $acc."_".$refbase."_sorted";
    my $sbam  = $acc."_".$refbase."_sorted.bam";

    my $cmd = "$bwa mem -t $threads $ref $f $r > $sam";
    run_cmd($cmd);
    undef $cmd;
    #say STDERR "===> Converting sam to bam and sorting...";

    ## sam -> bam
    $cmd = "$samtools view -bS $sam > $bam";
    run_cmd($cmd);
    undef $cmd;

    ## sort
    $cmd = "$samtools sort $bam $sort";
    run_cmd($cmd);

    return $sbam;
}

sub index_ref ($ref, $bwa, $samtools) {
    #say STDERR "===> Indexing reference $ref ...";
 
    ## index ref
    my $bcmd = "$bwa index $ref";
    my $scmd = "$samtools faidx $ref";
    
    run_cmd($bcmd);
    run_cmd($bcmd);
}


sub run_vep ($fvcf, $cvcf, $species) {
    #say STDERR "===> Running VEP ...";
    my $cvep_out  = $cvcf;
    my $fvep_out  = $fvcf;
    $cvep_out =~ s/\.vcf//;
    $fvep_out =~ s/\.vcf//;
    $cvep_out .= "_vep-out.txt";
    $fvep_out .= "_vep-out.txt";

    my $cvep_stats = $cvep_out."_summary.html";
    my $fvep_stats = $fvep_out."_summary.html";
    unlink $cvep_out if -e $cvep_out;
    unlink $fvep_out if -e $fvep_out;
    unlink $cvep_stats if -e $cvep_stats;
    unlink $fvep_stats if -e $fvep_stats;

    my $exe = "perl -I/home/statonse/apps/ensembl-tools-release-78/scripts/variant_effect_predictor";
    my $vep = "/home/statonse/apps/ensembl-tools-release-78/scripts/variant_effect_predictor/";
    $vep .= "variant_effect_predictor.pl";
    
    my $cmd = "$exe $vep --offline -i $fvcf --species $species -o $fvep_out";
    run_cmd($cmd);
    undef $cmd;
    
    $cmd = "$exe $vep --offline -i $cvcf --species $species -o $cvep_out";
    run_cmd($cmd);
    #say STDERR "===> Done.";
}

sub run_cmd ($cmd) {
    my ($stdout, $stderr, @res) = capture { system([0..5], $cmd); };

    #my @job
    #try {
	#@job = capture([0..5], $cmd);
    #}
    #catch {
	#say "\nERROR: $cmd exited. Here is the exception: $_\n";
	#exit;
    #};
}

