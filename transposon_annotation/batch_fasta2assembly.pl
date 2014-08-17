#!/usr/bin/env perl

# TODO: 

use strict;
use warnings;
use Getopt::Long; 
use File::Copy;
use IO::Compress::Bzip2 qw(bzip2 $Bzip2Error);

my $usage = "batch_fasta2assembly.pl -i /path/to/dir/of/fasta/files <options>

\tThe files in the input directory must end with .fasta for this
\tscript to work.

\tOptions: 
\t --clean     :    Remove all assembly files, keeping only the *Metrics.txt and *LargeContigs.fna
\t                  in a separate directory (renamed).
\t --compress  :    Compress the output directories.
\t --large     :    Return only the large contigs (>500 bp).
\t --all       :    Return all of the contigs.
\t -q|quiet    :    Don't print what Newbler is doing with the assemblies.
\t -v|verbose  :    Print assembly progress and messages to screeen.\n";

my $in_dir;
my $clean;
my $compress;
my $large;
my $all;
my $quiet;
my $verbose;

GetOptions(
	   "i|indir=s"    => \$in_dir,
	   "clean"        => \$clean,
	   "compress"     => \$compress,
	   "large"        => \$large,
	   "all"          => \$all,
	   "q|quiet"      => \$quiet,
	   "v|verbose"    => \$verbose,
	    );

if ($in_dir) {
    if ($verbose) {
	print "\n============= Starting assemblies ...\n";
    }
} 
else {
    die "\nERROR: No input directory was given!\n\n",$usage,"\n";
} 

if (!$large && !$all) {
    die "\nERROR: The options --large or --all must be specified!\n\n",$usage,"\n";
}

# from batch_genscan
my $date = `date`;
chomp $date;
$date =~ s/\s+/\_/g;
$date =~ s/\:/\-/g;     # archiving the directory is not possible with ":" as a delimiter in the date

my $assemb_out = "cluster_assemblies-".$date;
my $Metrics = "cluster_assembly_metrics-".$date;

opendir( DIR, $in_dir ) || die "Can't open directory:\n$in_dir\n"; 
my @fasta_files = grep /\.fasta$/, readdir DIR ; # using only .fasta right now
closedir( DIR );
my $num_files = @fasta_files;

# Start the program or die with a usage statement if no infiles found
if ($num_files > 0) {
    if ($verbose) {
	print "\n\n", "============= Running batch_fasta2assemb ============= ", "\n\n";
    }
}
else {
    die "\n","ERROR: No infiles could be found\n\n",$usage,"\n\n"; 
}

chdir($in_dir);

unless (-e $assemb_out) {    
    mkdir($assemb_out) || die "Could not create the output directory:\n$assemb_out\n";
}
 
unless (-e $Metrics) {
    mkdir($Metrics) || die "Could not create the output directory:\n$Metrics\n";
}

for my $ind_file (@fasta_files) {
    my $infile = $ind_file;
    $ind_file =~ s/\.fasta$//;

    my $proj_dir = $ind_file ."_newbler_out";
    if ($quiet) {
	my $assemb_cmd = "runAssembly -o $proj_dir -cpu 2 $infile 2>&1 > /dev/null";
	system($assemb_cmd); #print "\n";
    } 
    else {
	my $assemb_cmd = "runAssembly -o $proj_dir -cpu 2 $infile";
	system($assemb_cmd); print "\n";
    }
    chdir($proj_dir);
    if ($large) {
	if (-e "454LargeContigs.fna") {
    
	    my $largectgs = "454LargeContigs.fna";
	    my $assem = "454LargeContigs.fna";
	    my $metrics_file = "454NewblerMetrics.txt";
	    my $assem_stats = "454NewblerMetrics.txt";
	    $assem =~ s/454LargeContigs.fna/$ind_file/;
	    $assem .= "_" . "454LargeContigs.fna";
	    $assem_stats =~ s/454NewblerMetrics.txt/$ind_file/;
	    $assem_stats .= "_" . "454NewblerMetrics.txt";

	    copy("$largectgs","../$assemb_out/$assem") || die "Copy failed: $!";
	    copy("$metrics_file","../$Metrics/$assem_stats") || die "Copy failed: $!";

	}
    }
    if ($all) {
	if (-e "454AllContigs.fna") {
    
	    my $allctgs = "454AllContigs.fna";
	    my $assem = "454AllContigs.fna";
	    my $metrics_file = "454NewblerMetrics.txt";
	    my $assem_stats = "454NewblerMetrics.txt";
	    $assem =~ s/454AllContigs.fna/$ind_file/;
	    $assem .= "_" . "454AllContigs.fna";
	    $assem_stats =~ s/454NewblerMetrics.txt/$ind_file/;
	    $assem_stats .= "_" . "454NewblerMetrics.txt";

	    copy("$allctgs","../$assemb_out/$assem") || die "Copy failed: $!";
	    copy("$metrics_file","../$Metrics/$assem_stats") || die "Copy failed: $!";

	}
    }
	
    #print "\n$assem_stats\n","\n$assem\n";
    #print "\n$assemb_out\n","\n$Metrics\n";
    
    #copy("$largectgs","../$assemb_out/$assem") || die "Copy failed: $!";
    #system("cp $proj_dir/$largectgs ../../$assemb_out/$assem");
    #copy("$metrics_file","../$Metrics/$assem_stats") || die "Copy failed: $!";
    #system("cp $proj_dir/$metrics_file ../../$Metrics/$assem_stats");
    

    chdir("..");    
    # now clean up, could make an option
    if ($clean) {
	system("rm -rf $proj_dir/sff");
	system("rm -rf $proj_dir");
    }    

    #unlink($infile);     # not using a copy now!
    #chdir("..");
    if ($compress) {
	#my $cmpr_cmd = "tar -cjf $proj_dir.tar.bz2 $proj_dir";
	#system($cmpr_cmd);
	#system("rm -rf $proj_dir");
	bzip2 '<$proj_dir/454*>' => '<*bzip2>' || die "bzip2 failed: $Bzip2Error\n";
    }
}

if ($verbose) {
    print "\n\n", "============= batch_fasta2assemb Finished ============= ", "\n\n";
}

