#!/usr/bin/perl -w
#_______________________________________________________________________+
#                                                                       |
# batch_sixpack.pl               
#_______________________________________________________________________+
#                                                                       |
# Description: 
#                                                                       |
# Author: S. Evan Staton                                                |
# Contact: statonse<at>uga.edu                                          |
# Started: 3.3.11                                                      |                                                
# Updated:                                                      
#                                                                       |
# Suggested usage:                                                      |

#_______________________________________________________________________+
# TODO: 


use strict;
use Getopt::Long;
use File::Basename;
use File::Copy;

my $indir; # = $ARGV[0] || die "\nERROR: No input directory found!\n",$usage;
my $outdir;
my $orflen;
my $clean;
my $len;

# sequence counters
my $seqstot = 0;
my $transseqstot = 0;

# base counters
my $a = 0;
my $t = 0;
my $g = 0;
my $c = 0;
my $e = 0; 

GetOptions(
	   "i|indir=s"            => \$indir,
	   "o|outdir=s"           => \$outdir,
	   "orflen|orflength=i"   => \$orflen,
	   "clean"                => \$clean,
	  );
  
if (!$indir && !$outdir) {
    print "\nERROR: No input was given.\n";
    &usage();
    exit(0);
}

opendir( DIR , $indir ) || die "\nERROR: Could not open directory: $indir\n";
my @fastas = grep /\.fasta$|\.fas$|\.fna$/, readdir DIR;
closedir(DIR);

my $fasnum = @fastas;
if ($fasnum >= 1) {
    print "\n========== Translating sequences with minimum ORF length of $orflen.\n";
} else {
    die "\nERROR: No fasta files were found! Must end with .fasta, .fas, or .fna.\n";
}

chdir($indir);

foreach my $fas (@fastas) {
    $seqstot++;
    #my ($file,$dir,$ext) = fileparse($fas, qr/\.[^.]*/);   # parsing won't work if there are delimiters in the header
    
    # simple fileparse
    my $file = $fas;
    $file =~ s/\.fa.*$//;
    if ($file =~ m/\:/g) {
	warn "\nERROR: $file contains delimiters (i.e. : ) in the header for $fas. Exiting.\n";
	next;
    } elsif ($file =~ m/\;/g) {
	warn "\nERROR: $file contains delimiters (i.e. ; ) in the header for $fas. Exiting.\n";
	next;
    } else {

	open( FAS , $fas ) || die "\nERROR: Could not open file: $fas\n";
	while(<FAS>) {
	    chomp;
	    unless(/^\>/) {
		my @seq;
		push(@seq,$_);
		my $dna = join(' ',@seq);
		$dna =~ s/\s//g;
		$len = length($dna);
		while($dna =~ /a/ig){$a++;}
		while($dna =~ /t/ig){$t++;}
		while($dna =~ /g/ig){$g++;}
		while($dna =~ /c/ig){$c++;}
		while($dna =~ /[^atgc]/ig){$e++;} # could store non-atgc bases in a hash, sort and count them
                                                  # but that's not really necessary and will just make this slower
		unlink($dna);
	    } 
	}
	close(FAS);
	if ($e > 0 ) {
	    warn "\nWARNING: $fas (length: $len) contains $e non-nucleic symbols. Skipping.\n";
	    next;
	    unlink($e);
	} else {
	    unlink($a,$t,$g,$c);
	    my $pepfile = $file.".faa";
	    my $transout = $file.".sixpack";
	    my $transcmd = "sixpack ".
		           "-sequence $fas ".
	                   "-outseq $pepfile ".
		           "-outfile $transout ".
		           "--orfminsize $orflen ". 
		           "-debug ".
			   #"-auto ".     # be quiet unless ...
			   "-error ".    # there are errors, then report them and ...
			   "-die ";      # any dying messages
                          # will want this to run quietly when I get all the formatting fixed
	    system($transcmd);
	    if (-s $pepfile) {
		$transseqstot++;
		move("$pepfile","../$outdir") || die "Copy failed: $!";
	    }
	    if ($clean) {
		unlink($transout);
	    }
	}
    }
}

print "\n========== $seqstot sequences processed.\n";
print "\n========== $transseqstot sequences translated with ORFs above $orflen.\n";

exit;

#---------------------
# Subroutines
#---------------------

sub usage {
    my $script = basename($0);
    print STDERR <<END
USAGE: $script -i indir -o outdir 

\tRequired:
\t -i|indir            :       The input directory of nucleotide fasta* files.
\t                             *Must end with .fasta, .fas, or .fna at this time.
\t -o|outdir           :       A place to put the translated sequences.

\tOptions:
\t -orflen|orflength   :       An interger that will serve as the lower threshold
\t                             length for ORFs to consider prior to translating.
\t --clean             :       Remove all of the alignment files created by EMBOSS
\t                             if only the translated sequences are desired.\n
END
}
