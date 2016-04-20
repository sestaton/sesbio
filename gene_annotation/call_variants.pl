#!/usr/bin/env perl

use 5.020;
use strict;
use warnings;
use autodie;
use File::Spec;
use File::Find;
use File::Basename;
use Parallel::ForkManager;
use Try::Tiny;
use IPC::System::Simple qw(system);
use Capture::Tiny       qw(:all);
use List::MoreUtils     qw(natatime);
use experimental 'signatures';

my @reads;
my $usage    = "$0 ref\n";
my $ref      = shift or die $usage;
my $species  = 'scf7180038271797';
my $dir      = 'cleaned_transcriptome_reads';
my $github   = File::Spec->catfile($ENV{HOME}, 'github');
my $bwa      = File::Spec->catfile($github, 'bwa', 'bwa');
my $samtools = File::Spec->catfile($github, 'samtools', 'samtools');
my $bcftools = File::Spec->catfile($github, 'bcftools', 'bcftools');
my $vcfutils = File::Spec->catfile($github, 'bcftools', 'vcfutils.pl');
my $threads  = 12;
my $excl     = 'EPSPS_reads|unpe|test';
#index_ref($ref, $bwa, $samtools);

find( sub { push @reads, $File::Find::name if -f and /\.fq$/ }, $dir );
my @rawreads = grep { /[12]/ } 
               grep { ! /$excl/ } @reads;

my @sreads = sort @rawreads;
my $it = natatime 2, @sreads;

my $pm = Parallel::ForkManager->new($threads);
$pm->run_on_finish( sub { my ($pid, $exit_code, $ident, $exit_signal, $core_dump) = @_;#
			  say basename($ident)," just finished with PID $pid and exit code: $exit_code";
		    } );

while (my @vals = $it->()) {
    my ($f, $r) = @vals;
    my ($acc) = ($f =~ /[weedy|wild|elite]\.(.*)\.[12]\.fq$/);
    $pm->start($acc) and next;

    my $bam = map_reads($ref, $acc, $bwa, $f, $r, $threads);
    my ($fvcf, $cvcf) = call_variants($ref, $bam, $samtools, $bcftools, $vcfutils);
    run_vep($fvcf, $cvcf, $species);

    $pm->finish;
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

    ##NB: this is kind of a quick-and-dirty method, replace with GATK
    ##    RNA-Seq pipeline for more important analyses
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

    my $veplib = File::Spec->catdir($ENV{HOME},'apps','ensembl-tools-release-78','scripts','variant_effect_predictor');
    my $exe = "perl -I$veplib";
    my $vep = File::Spec->catfile($veplib, 'variant_effect_predictor.pl');
    
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

