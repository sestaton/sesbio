#!/usr/bin/env perl

=head1 NAME 
                                                                       
 gsAssembler.pl - Sample a flowgram n times and do n assemblies          

=head1 SYNOPSIS    

 gsAssembler.pl -i reads.sff -p job_name -o assem_stats.txt

=head1 DESCRIPTION
                                                                   
 Takes as input a flowgram file and does a de novo        
 assembly on, optionally, a subset of the reads with the 
 program gsAssembly. The number of samples taken is an option
 specified at the command line. Will optionally take n samples 
 of i bases or j reads.    
                                                                       
 This could be useful for fine tuning the assembly parameters          
 by keeping results of each run in a spreadsheet or for general        
 record keeping.                                                        
                                                                       
=head1 AUTHOR 

S. Evan Staton                                                

=head1 CONTACT
 
evan at evanstaton dot com

=head1 REQUIRED ARGUMENTS

=over 2

=item -i, --infile

The flowgram to be used for assembly.

=item -p, --prefix

A short prefix to be used to name the assembly output directory.

=back

=head1 OPTIONS

=over 2

=item -t, --trimfile

A file to use for contaminant filtering.

=item -o, --outfile

A file to hold the assembly statistics. 

=item --sample

This option is required if a subsample of the flowgram is to be generated.
Required if --bases or --reads is specified.

=item --bases

Allows one to specify the number of bases to be taken in the form 1m for
1 million or 1k for 1 thousand. 

=item --reads

Allows one to specify the number of reads to be taken in the form 5000k for
5 thousand reads. 

=item -s, --sample_cov

The desired coverage of the genome or BAC.

=item -n,--num_subsamples

The number of samples to generate at the specified level.

=item --consed

Generate the directories and files for launching Consed after assembly.

=cut                                          

#
# TODO: Add example usage to POD
#       Add try/catch, ipcss, and method for handling commands
#
use 5.010;
use strict;
use warnings;
use File::Basename;
use File::Copy;
use Getopt::Long;
use Pod::Usage;

my $infile;
my $outfile;
my $prefix;
my $sample;
my $bases;
my $reads;
my $sample_cov;
my $num_subsamples;
my $consed;
my $trimfile;

GetOptions(#Required arguments 
	   'i|infile=s'         => \$infile,
	   'o|outfile=s'        => \$outfile,
	   #Options
	   't|trimfile=s'       => \$trimfile,
	   'p|prefix=s' 	=> \$prefix,
	   'sample'             => \$sample,
	   'bases'      	=> \$bases,
	   'reads'      	=> \$reads,
	   's|sample_cov=s'     => \$sample_cov,
	   'n|num_subsamples=i' => \$num_subsamples,
	   'consed'             => \$consed,
	  );

# Start the program or die with a usage statement and specific error message
# check for infile
if (!$infile || !$prefix) {
    say "\nERROR: No input file or prefix were given.";
    usage();
    exit(0);
}
else {
    say "\n============= Generating subsamples for file: ",
    $infile if $sample;
    say "\n============= Generating assembly for file: ",
    $infile if !$sample;
}

#check for prefix
if ($prefix) {
    say "============= Generating a subsample of reads for project name: ",
    $prefix if $sample;
    say "============= Generating assembly for project name: ",
    $prefix if !$sample;
}

#check for sample_cov
if ($sample) {
    if ($sample_cov) {
	say "============= Generating a subsample of reads for projected coverage: ",
	$sample_cov;
    }
    else {
	say "\nERROR: A sample coverage was not given at the command line.";
	usage();
	exit(0);
    }
    
    #check for num_subsamples
    if ($num_subsamples) {
	say "============= Generating ", $num_subsamples,
	" subsamples of reads for file: ", $infile;
    }
    else {
	say "\nERROR: The number of subsamples was not given at the command line.";
	usage();
	exit(0);
    }
    
    #check for --bases or --reads
    if ( $bases && $reads ) {
	say "\nERROR: The --reads or --bases options must not be given together at the command line.";
	usage();
	exit(0);
    }
    elsif ( !$bases && !$reads ) {
	say "\nERROR: The --reads or --bases options must be given at the command line.";
	usage();
	exit(0);
	}
    elsif ($bases) {
	say "============= Generating subsamples of bases for: ", $sample_cov,
	" bases\n";
    }
    elsif ($reads) {
	say "============= Generating subsamples of reads for: ", $sample_cov,
	" reads\n";
    }
}

#-----------------------------------------------------+
# Take an input sff and generate a subsample of reads |
#-----------------------------------------------------+
#find the sfffile program
my $sfffile = `which sfffile`;
chomp $sfffile;

#system("$sfffile") == 0 || die "Could find the program sfffile: Please make sure it is installed and in your $PATH";

my $results_dir = $prefix . "_newbler_split";
mkdir $results_dir
  || die "Could not create directory to write results\n";   #chdir $results_dir;

if ($sample) {
    chdir $results_dir ;
    my $it = 0;    # initialize count of subsample iterations
    for ( $it = 1; $it <= $num_subsamples; $it++ ) { # this is not very Perlish, use range instead
	
	my $its = ($it);
	say "\n============> Generating subsample ", $its,
	" of reads for file: $infile\n";

	if ($reads) {
	    my $sff_cmd_reads = "sfffile -pickr $sample_cov -o sub_$its-$infile ../$infile";
	    system($sff_cmd_reads) == 0 or die $!;
	    print "\n";
	}	
	elsif ($bases) {
	    my $sff_cmd_bases = "sfffile -pickb $sample_cov -o sub_$its-$infile ../$infile";    
	    # 4m = 32x for a 125kb BAC
	    system($sff_cmd_bases) == 0 or $!;
	    print "\n";
	}
	if ( $it == $num_subsamples ) {
	    say "\n\n", "============= Done sampling reads in: ", $infile,
	    " for ", $num_subsamples, " subsamples =============\n";
	}

    }
    #system("mv sub_* $results_dir");    # safer to `mv stuff` ?
    #move("sub_*","$results_dir") || die "Move failed: $!";
}
#--------------------------------------+
# DO runAssembly on generated sfffiles |
#--------------------------------------+
# find runAssembly
my $runAssembly = `which runAssembly`;
chomp $runAssembly;

#system("$runAssembly") == 0 || "Couldn't find runAssembly: Please make sure it is installed and in your $PATH";

# get the sfffiles for assembly
if (!$sample) {
    chdir $results_dir;
    copy("../$infile",".") || die "Copy failed: $!";
}
my $results_path = `pwd`;
chomp $results_path;
$results_path =~ s/ \/[^\/]+$//;    # removes the last / and everything after it
chdir("..");
opendir my $dir, $results_dir || die "Can't open directory: $results_dir";
my @sub_sffs =
  map { "$results_path/$_" } grep { $_ !~ /^\./ } readdir($dir);
closedir $dir;

# do the assembly on each subsample
for my $f (@sub_sffs) {
    unless ( ( $f eq "." ) || ( $f eq ".." ) ) {
	say "\n============> Running Newbler on file: $f\n";

	# cleaned up the output directory names 11/02/10 SES
	my $output = $f;
	$output =~ s/\.[^\.]*$//;    # http://www.perlmonks.org/?node_id=729477
	if ($consed) {
	    my $run_Assemb_consed_cmd;
	    $run_Assemb_consed = "$runAssembly -o $output -consed -vt $trimfile $f" if $trimfile;
	    $run_Assemb_consed = "$runAssembly -o $output -consed $f" if !$trimfile;
	    system($run_Assemb_consed) == 0 or die $!;
	    print "\n";
	}
	else {
	    my $run_Assemb = "$runAssembly -o $output -vt $trimfile $f";
	    system($run_Assemb) == 0 or die $!;
	    print "\n";
	}
	if ($outfile) {
	    chdir $output;
	    #my $sffsample = $output.".sff";
	    #my $outparsed = $output."_parsed.txt";
	    print "\n============= Parsing assembly output for $f\n";
	    my $Metrics = "454NewblerMetrics.txt";
	    get_NewblerMetrics($Metrics,$outfile,$output);
	    chdir("..");
	}
    }
}

system("cat *_parsed.txt > $outfile"); # shame! re-write this in perl
unlink $infile;
system("rm *_parsed.txt");             # more shame

if ($sample){
    say "\n============= Assembly complete for each of ", $num_subsamples,
    " subsamples for file: ", $infile, " =============\n";
} 
else {
    say "\n============= Assembly complete for file: $infile  =============\n";
}
#-------------
# SUBROUTINES
#-------------
sub get_NewblerMetrics {
	my ( $infile, $outfile, $parsed ) = @_;
	$parsed .= "_parsed.txt";
	open my $in, '<',  $infile or die "\nERROR: Could not open file: $infile\n";
	open my $out,'>', $outfile or die "\nERROR: Could not open file: $outfile\n";

	say $out join "\t", "#filename","num_reads","num_bases","tot_num_reads","tot_num_bases",
			    "num_aln_reads","num_aln_bases","infer_read_error","num_assemb",
			    "num_partial","num_singleton","num_repeat","num_outlier","num_too_short",
			    "num_contigs","num_bases_large","ave_contig_size","N50_contig_size",
			    "largest_contig_size","Q40PlusBases","Q39MinusBases","num_contigs",
			    "num_bases_all";

	my $line = 0;
	while (<$in>) {
	    chomp;
	    $line++;
	    if (/(path = )(\"\/.*\";$)/ && $line == 20) {  # works now, matching only the first "path" 3/19/11 SES
		#if (/path = / && $line == 20) {    
		
		my $flowgram = $2;
		#my ($offset,$pathdir,$eq,$path) = split(/\s+/,$flowgram);
		#$path =~ s/\"//g;
		#$path =~ s/\;$//;
		$flowgram =~ s/\"//g;
		$flowgram =~ s/\;$//;
		my ( $sfffile, $dir, $ext ) = fileparse( $flowgram, qr/\.[^.]*/ );
		$sfffile .= $ext;
		print $out $sfffile."\t";
	    }
	    elsif (m/(numberOfReads = )(\d+)\, (\d+)\;/) {
		#numberOfReads = 10937, 10937;
		my $untrimmedReads = $2;
		my $trimmedReads   = $3;
		print $out $untrimmedReads."\t".$trimmedReads."\t";
	    }
	    elsif (m/(numberOfBases = )(\d+)\, (\d+)\;/) {
		#numberOfBases = 4000256, 3993688;
		my $untrimmedBases = $2;
		my $trimmedBases   = $3;
		print $out $untrimmedBases."\t".$trimmedBases."\t";
	    }
	    elsif (m/(totalNumberOfReads = )(\d+)\;/) {
		#totalNumberOfReads = 10937;
		my $totNumReads = $2;
		print $out $totNumReads."\t";
	    }
	    elsif (m/(totalNumberOfBases = )(\d+)\;/) {
		#totalNumberOfBases = 3993688;
		my $totNumBases = $2;
		print $out $totNumBases."\t";
		
		# Omitted for now
		#numberSearches   = 1534;
		#seedHitsFound    = 442889, 288.72;
		#overlapsFound    = 46357, 30.22, 10.47%;
		#overlapsReported = 42099, 27.44, 90.81%;
		#overlapsUsed     = 20442, 13.33, 48.56%;
	    }
	    elsif (m/(numAlignedReads     = )(\d+)\, (\d\d\.\d\d\%)\;/) {
		#numAlignedReads     = 10677, 97.62%;
		my $numAlnReads  = $2;
		my $percAlnReads = $3;
		print $out $numAlnReads."\t".$percAlnReads."\t";
	    }
	    elsif (m/(numAlignedBases     = )(\d+)\, (\d\d\.\d\d\%)\;/) {
		#numAlignedBases     = 3940332, 98.66%;
		my $numAlnBases  = $2;
		my $percAlnBases = $3;
		print $out $numAlnBases."\t".$percAlnBases."\t";
	    }
	    elsif (m/(inferredReadError  = )(\d\.\d\d\%)\, (\d+)/) {
		#inferredReadError  = 0.84%, 33107;
		my $percAlnError = $2;
		my $numAlnError  = $3;
		print $out $percAlnError."\t".$numAlnError."\t";
	    }
	    elsif (m/(numberAssembled = )(\d+)\;/) {
		#numberAssembled = 10234;
		my $numReadsAssem = $2;
		print $out $numReadsAssem."\t";
	    }
	    elsif (m/(numberPartial = )(\d+)\;/) {
		#numberPartial   = 443;
		my $numPartial = $2;
		print $out $numPartial."\t";
	    }
	    elsif (m/(numberSingleton = )(\d+)\;/) {
		#numberSingleton = 85;
		my $numSingleton = $2;
		print $out $numSingleton."\t";
	    }
	    elsif (m/(numberRepeat    = )(\d+)\;/) {
		#numberRepeat    = 7;
		my $numRepeat = $2;
		print $out $numRepeat."\t";
	    }
	    elsif (m/(numberOutlier   = )(\d+)\;/) {
		#numberOutlier   = 27;
		my $numOutlier = $2;
		print $out $numOutlier."\t";
	    }
	    elsif (m/(numberTooShort  = )(\d+)\;/) {
		#numberTooShort  = 141;
		my $numTooShort = $2;
		print $out $numTooShort."\t";
	    }
	    elsif (m/(numberOfContigs   = )(\d+)\;/) {
		#	numberOfContigs   = 8;
		my $numContigs = $2;
		print $out $numContigs."\t";
	    }
	    elsif (m/(numberOfBases     = )(\d+)\;/) {
		#	numberOfBases     = 109378;
		my $numBases = $2;
		print $out $numBases."\t";
	    }
	    elsif (m/(avgContigSize     = )(\d+)\;/) {
		#	avgContigSize     = 13672;
		my $aveContigSize = $2;
		print $out $aveContigSize."\t";
	    }
	    elsif (m/(N50ContigSize     = )(\d+)\;/) {
		#	N50ContigSize     = 35275;
		my $N50ContigSize = $2;
		print $out $N50ContigSize."\t";
	    }
	    elsif (m/(largestContigSize = )(\d+)\;/) {
		#	largestContigSize = 38319;
		my $largestContigSize = $2;
		print $out $largestContigSize."\t";
	    }
	    elsif (m/(Q40PlusBases      = )(\d+)\, (\d\d\.\d\d\%)\;/) {
		#	Q40PlusBases      = 109302, 99.93%;
		my $Q40PlusBaseCount = $2;
		my $Q40PlusBasePerc  = $3;
		print $out $Q40PlusBaseCount."\t".$Q40PlusBasePerc."\t";
	    }
	    elsif (m/(Q39MinusBases     = )(\d+)\, (\d\.\d\d\%)\;/) {
		#	Q39MinusBases     = 76, 0.07%;
		my $Q39MinusBaseCount = $2;
		my $Q39MinusBasePerc  = $3;
		print $out $Q39MinusBaseCount."\t".$Q39MinusBasePerc."\t"; 
	    }
	    elsif (m/(numberOfContigs = )(\d+)\;/) {
		#	numberOfContigs = 8;
		my $numContigs = $2;
		print $out $numContigs."\t";
	    }
	    elsif (m/(numberOfBases   = )(\d+)\;/) {
		#	numberOfBases   = 109378;
		my $numBases = $2;
		print $out $numBases."\n";
	    }
	}
	close $in;
	close $out;
	copy("$outfile", "$parsed") || die "Copy failed: $!";
	move("$parsed", "../")  || die "Move failed: $!";
}

sub usage {
	my $script = basename($0);
	print STDERR <<END
USAGE: $script [-i][-o][-n][-s][--sample][--reads][--bases][--consed]

USAGE examples:

    Ex.1 : Use all reads and generate assembly report.
	
    	  $script -i myreads.sff -p MID2 -o assemb_stats.txt 
           
    Ex.2 : Generate 3 samples of 4 million bases.

	  $script -i myreads.sff -p MID2 -o assem_stats.txt --sample \ 
	    		-s 4m -n 3 --bases --consed        
                          
    Ex.3 : Generate 3 samples from 7500 reads.                  
                                                                       
          $script -i myreads.sff -p MID2 -o assemb_stats.txt --sample \ 
        		-s 7500k -n 3 --reads --consed   
	 		
Required:
   -i|infile         :       The Roche flowgram.
                             
   -p|prefix         :       The name of the assembly job.

Options:
   -o|outfile        :       A file to hold the parsed output.
   -n|num_subsamples :       Number of subsamples to generate.
   -s|sample_cov     :       The target coverage.
   --sample          :       Sample the flowgram instead of 
                             using all reads. Must be used in
                             combination with <--reads> or
                             <--bases> to specify how the 
                             flowgram is to be sampled.
   --reads           :       Generate a sample of the specified 
                             number of reads in the form n(m|k)
  			     where n is an integer, m=millions,
  			     k=thousands.
   --bases           :	     Generate a sample of the specified 
                             number of bases in the form n(m|k)
  			     where n is an integer, m=millions,
  			     k=thousands.
   --consed          :       Create the directory structure and
                             files needed by Consed.\n                            
END
}

#-------------------+
# REVISION HISTORY  |
#-------------------+
# 5.15.10
# -Started the program to manage taking subsamples
# of reads from an input sff and running newbler
# on each. Right now just takes a subsample of
# reads
# 5.16.10
# -program now takes the requested number of samples
# and assembles each sample. There is a pesky
# warning dealing with how the files are being named:
# "Useless use of a constant in void context" with how $it
# is being used to get sfffile counts
# Otherwise, all is working
# 5.17.10
# -rewrote sampling routine so the the constant being
# used is in context. Changed the name of the program
# to sub_runRA.pl
# added more options and started working on getting
# assembly statistics. Error handling needs attention
# because program does not exit cleanly if some
# required options are missing. DONE 3:14 am
# 3.15.11
# - added parsing subroutine, POD, option to just
# run the assembly and parse the output without
# generating samples
# 8.22.13
# - clean up code, use modern pragmas

# Started: 5.15.10                                                      |
# Updated: 5.17.10                                                      |
#                                                                       |
# Suggested usage:                                                      |
#   sample_n_runAssemb.pl -i infile -p job_name -s sample_coverage      |
#                         -n num_subsamples --bases n <options>         |
# Examples:                                                             |
#   sub_runRA.pl -i myreads.sff -p MID2 -s 4m -n 3                      |
#                          --bases 4m --consed                          |
#                                                                       |
#   sub_runRA.pl -i myreads.sff -p MID2 -s 4m -n 3                      |
#                          --reads 7500k --consed                       |
#                                                                       |
# Notes: The 4m above means 4 million bases and 7500k means 7500 reads  |
#        m=millions, k=thousands                                        |
#        Either --bases or --reads must be given but not together       |
#                                                                       |
# Requires: Newbler software installed and in your $PATH                |
#                                                                       |
#_______________________________________________________________________+
